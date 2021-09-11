import sys
from functools import *
import trimesh
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
    file.write("comment https://github.com/mikedh/trimesh\n")
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

def produce_slices(mesh, _dir='x', n_slices=20, out_filename = 'output/slice'):
  bounds = trimesh.bounds.corners(mesh.bounding_box.bounds)
  min_corner = bounds[0,:]
  max_corner = bounds[6,:]

  direction = np.array([0, 0, 0])
  if _dir == 'x':
    direction[0] = 1
    step_size = (max_corner[0] - min_corner[0]) / (n_slices + 1)
  elif _dir == 'y':
    direction[1] = 1
    step_size = (max_corner[1] - min_corner[1]) / (n_slices + 1)
  else:
    direction[2] = 1
    step_size = (max_corner[2] - min_corner[2]) / (n_slices + 1)

  start = min_corner + (0.5 * step_size) * direction
  for n in range(n_slices):
    origin = start + (n * step_size) * direction
    path3d = slice_mesh(mesh, direction, origin)
    export_polyline(path3d, "polyline_" + str(n) + ".ply")
    path, _ = path3d.to_planar()
    filename = out_filename + '_' + str(n) + '.svg'
    trimesh.path.exchange.export.export_path(path, file_type='svg', file_obj=filename)

if __name__ == "__main__":
  produce_slices(trimesh.load_mesh(sys.argv[1]))
