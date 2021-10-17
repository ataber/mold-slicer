import sys
import os
import argparse
from math import pi
from loguru import logger
from functools import *
import trimesh
import trimesh.transformations as transformations
import numpy as np
import random
from shapely.geometry import LineString

def slice_mesh(mesh, normal, origin):
  return mesh.section(plane_origin=origin, plane_normal=normal)

def export_polyline(path3d):
  polylines = []
  for entity in path3d.entities:
    polyline = []
    for index in entity.points:
      polyline.append(path3d.vertices[index])
    polylines.append(polyline)
  return polylines

def write_polyline_ply(polylines, filename):
  n_vertices = reduce(lambda a,b: a+b, [(len(line) - 1) * 3 for line in polylines])
  n_faces = reduce(lambda a,b: a+b, [(len(line) - 1) for line in polylines])
  with open(filename, "w") as file:
    file.write("ply\n")
    file.write("format ascii 1.0\n")
    file.write("element vertex " + str(n_vertices) + '\n')
    file.write("property float x\n")
    file.write("property float y\n")
    file.write("property float z\n")
    file.write("element face " + str(n_faces) + '\n')
    file.write("property list uchar int vertex_indices\n")
    file.write("end_header\n")
    vertices = []
    faces = []
    for polyline in polylines:
      for v in range(len(polyline) - 1):
        i1 = len(vertices)
        vertices.append(polyline[v])
        i2 = len(vertices)
        vertices.append(polyline[v + 1])
        i3 = len(vertices)
        vertices.append(polyline[v])
        faces.append([i1, i2, i3])
    for vertex in vertices:
      file.write("{} {} {}\n".format(*vertex))
    for face in faces:
      file.write("3 {} {} {}\n".format(*face))

def rotation_matrix_from_vectors(vec1, vec2):
  """ Find the rotation matrix that aligns vec1 to vec2
  :param vec1: A 3d "source" vector
  :param vec2: A 3d "destination" vector
  :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
  """
  a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
  v = np.cross(a, b)
  c = np.dot(a, b)
  s = np.linalg.norm(v)
  kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
  rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
  return rotation_matrix

def produce_slices(mesh, step_size=0.1875, out_dir='output', negative=False):
  logger.warning(out_dir)
  # first clear the output directory
  os.system(f"""rm -rf {out_dir} && mkdir {out_dir} && mkdir {out_dir}/polyline && mkdir {out_dir}/extrusion && mkdir {out_dir}/slice""")
  direction = np.array([0, 0, 0])
  direction[2] = 1

  bnds = trimesh.bounds.corners(mesh.bounding_box.bounds)
  min_corner = bnds[0,:]
  max_corner = bnds[6,:]

  if negative:
    transform = transformations.translation_matrix((max_corner + min_corner) / 2)
    box = trimesh.primitives.Box(extents=(max_corner - min_corner) * 1.05, transform=transform)
    # reverse orientation
    for face in mesh.faces:
      tmp = face[0]
      face[0] = face[1]
      face[1] = tmp
    mesh_to_slice = trimesh.util.concatenate([mesh, box])
  else:
    mesh_to_slice = mesh

  start = min_corner + (0.5 * step_size) * direction
  slice_plane_origin = start
  n = 0
  while True:
    slice_plane_origin = start + (n * step_size) * direction
    if slice_plane_origin[2] > max_corner[2]:
      return

    # export polyline
    path3d = slice_mesh(mesh_to_slice, direction, slice_plane_origin)
    polylines = export_polyline(path3d)
    write_polyline_ply(polylines, f'{out_dir}/polyline/polyline_{n}.ply')

    # export SVG
    path, _ = path3d.to_planar()
    slice_filename = f'{out_dir}/slice/slice_{n}.svg'
    trimesh.path.exchange.export.export_path(path, file_type='svg', file_obj=slice_filename)
    with open(slice_filename, 'r') as f:
      lines = f.read().splitlines()
    # set units
    with open(slice_filename, 'w') as f:
      for line in lines:
        if "width=" in line:
          f.write(f'width="{path.extents[0]}in"')
        elif "height=" in line:
          f.write(f'height="{path.extents[1]}in"')
        else:
          f.write(line)
        f.write('\n')

    # export extrusion
    centroid = np.zeros(3)
    for v in path3d.vertices:
      centroid += v
    centroid /= len(path3d.vertices)
    centroid[2] = min_corner[2]
    T = transformations.translation_matrix(centroid + (n * step_size) * direction)

    extrusion = path.extrude(step_size)
    if type(extrusion) is list:
      extrusion = trimesh.util.concatenate(extrusion)
    extrusion.apply_transform(T)
    trimesh.exchange.export.export_mesh(extrusion, file_type='ply', file_obj=f'{out_dir}/extrusion/extrusion_{n}.ply')
    n += 1

def rotation_matrix_from_vectors(vec1, vec2):
  a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
  v = np.cross(a, b)
  c = np.dot(a, b)
  s = np.linalg.norm(v)
  kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
  rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
  return rotation_matrix

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--files", help="list of files, comma separated")
  parser.add_argument("--up_dirs", help="list of up_dirs ('X', 'Y', or 'Z') per file, comma separated")
  args = parser.parse_args()
  for i, file in enumerate(args.files.split(',')):
    orientation = args.up_dirs.split(',')[i]
    name = file.split('.')[0]
    mesh = trimesh.load_mesh(file)
    if orientation == 'X':
      angles = [0, pi/2, 0]
    elif orientation == 'Y':
      angles = [pi/2, 0, 0]
    else:
      angles = [0, 0, 0]
    M = transformations.compose_matrix(translate=-mesh.centroid, angles=angles)
    mesh.apply_transform(M)
    produce_slices(mesh, out_dir=f'{name}_output')
