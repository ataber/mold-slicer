import sys
import trimesh
import numpy as np
import random
from shapely.geometry import LineString

def slice_mesh(mesh, normal, origin):
  return mesh.section(plane_origin=origin, plane_normal=normal)

def export_polyline(path3d, filename):
  polyline = []
  v = 0
  for entity in path3d.entities:
    for index in entity.points:
      polyline.append(path3d.vertices[index])
      v += 1

  write_polyline_ply(polyline, filename)

def write_polyline_ply(polyline, filename):
  with open(filename, "w") as file:
    file.write("ply")
    file.write("format ascii 1.0")
    file.write("comment https://github.com/mikedh/trimesh")
    file.write("element vertex " + str(len(polyline) * 3))
    file.write("property float x")
    file.write("property float y")
    file.write("property float z")
    file.write("element face " + str(len(polyline) - 1))
    file.write("property list uchar int vertex_indices")
    file.write("end_header")
    v = 0
    faces = []
    vertices = []
    while v < len(polyline) - 1:
      i1 = len(vertices)
      i2 = len(vertices)
      file.write("%d %d %d".format(*polyline[v + 1]))
      i3 = len(vertices)
      file.write("%d %d %d".format(*polyline[v]))
      faces.append([i1, i2, i3])
      v += 1
    for vertex in vertices:
      file.write("%d %d %d".format(*vertex))
    for face in faces:
      file.write("3 %i %i %i".format(*face))

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
    abirt
    # trimesh.path.exchange.export.export_path(path, file_type='svg', file_obj=filename)

if __name__ == "__main__":
  produce_slices(trimesh.load_mesh(sys.argv[1]))
