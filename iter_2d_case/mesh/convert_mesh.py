import meshio
mesh = meshio.read("Mesh_ITER.med")
#meshio.write(mesh_file_xdmf, mesh)  # won't work for FEniCS, mixed

# In order to use MeshFunction of FEniCS
# The tag must be a positive number (size_t)
mesh.cell_data_dict["cell_tags"]["triangle"] *= -1
mesh.cell_data_dict["cell_tags"]["line"] *= -1


print(mesh.cell_tags)
# mesh.cell_tags = {-6: ['Down'], -7: ['Top'], -8: ['Lying on Top']}


# Export mesh that contains only triangular faces
# along with tags
meshio.write_points_cells(
    "mesh_domains.xdmf",
    mesh.points,
    [mesh.cells[1]],
    cell_data={"f": [-1*mesh.cell_data["cell_tags"][1]]},
)

# Export mesh that contains only lines
# along with tags
meshio.write_points_cells(
    "mesh_boundaries.xdmf",
    mesh.points,
    [mesh.cells[0]],
    cell_data={"f": [-1*mesh.cell_data["cell_tags"][0]]},
)
