import sys
from functools import *
import trimesh
import trimesh.transformations as transformations
import numpy as np
import random
from shapely.geometry import LineString

def slice_mesh(mesh, normal, origin):
  return mesh.section(plane_origin=origin, plane_normal=normal)

def export_polyline(path3d, filename):
  polylines = []
  for entity in path3d.entities:
    polyline = []
    for index in entity.points:
      polyline.append(path3d.vertices[index])
    polylines.append(polyline)

  write_polyline_ply(polylines, filename)

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

def produce_slices(mesh, _dir=2, n_slices=40, out_dir='output', negative=True):
  bnds = trimesh.bounds.corners(mesh.bounding_box.bounds)
  min_corner = bnds[0,:]
  max_corner = bnds[6,:]

  direction = np.array([0, 0, 0])
  direction[_dir] = 1
  step_size = (max_corner[_dir] - min_corner[_dir]) / n_slices

  if negative:
    transform = transformations.translation_matrix((max_corner + min_corner) / 2)
    box = trimesh.primitives.Box(extents=(max_corner - min_corner), transform=transform)
    # reverse orientation
    for face in mesh.faces:
      tmp = face[0]
      face[0] = face[1]
      face[1] = tmp
    mesh_to_slice = trimesh.util.concatenate([mesh, box])
  else:
    mesh_to_slice = mesh

  start = min_corner + (0.5 * step_size) * direction
  for n in range(n_slices):
    slice_plane_origin = start + (n * step_size) * direction
    path3d = slice_mesh(mesh_to_slice, direction, slice_plane_origin)
    export_polyline(path3d, f'{out_dir}/polyline_{n}.ply')
    path, _ = path3d.to_planar()
    trimesh.path.exchange.export.export_path(path, file_type='svg', file_obj=f'{out_dir}/slice_{n}.svg')

    R = np.zeros((4, 4))
    if _dir != 2:
      rot = rotation_matrix_from_vectors(np.array([0, 0, 1]), direction)
      R[0:3,0:3] = rot
      R[3,3] = 1
    else:
      for i in range(4):
        R[i,i] = 1 # Identity

    centroid = np.zeros(3)
    for v in path3d.vertices:
      centroid += v
    centroid /= len(path3d.vertices)
    centroid[_dir] = min_corner[_dir]
    T = transformations.translation_matrix(centroid + (n * step_size) * direction)
    M = transformations.concatenate_matrices(R, T)

    extrusion = path.extrude(step_size)
    if type(extrusion) is list:
      extrusion = trimesh.util.concatenate(extrusion)
    extrusion.apply_transform(M)
    trimesh.exchange.export.export_mesh(extrusion, file_type='ply', file_obj=f'{out_dir}/extrusion_{n}.ply')

if __name__ == "__main__":
  produce_slices(trimesh.load_mesh(sys.argv[1]))
