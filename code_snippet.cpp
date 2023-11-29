 
int n_cells = solid_3d.triangulation.n_active_cells();

Eigen::MatrixXd cell_contributions = Eigen::MatrixXd::Zero(n_cells, n_cells);
        //std::vector<bool> visited_cells(n_cells, false);
bool *visited_cells=new bool[n_cells]();
for (const auto &cell : solid_3d.dof_handler_ref.active_cell_iterators())
  {

        if(cell->at_boundary()){
          continue;
          }

          const unsigned int cell_index = cell->active_cell_index();
          if (visited_cells[cell_index])
             continue; // Skip if the cell has been visited
            // Check if the cell has been visited
            // Loop over all degrees of freedom on the current cell
        for (unsigned int i = 0; i < solid_3d.dofs_per_cell; ++i)
      {
        //dealii::DoFHandler::cell_iterator neighbor;
        // Loop over all neighboring cells
        for (unsigned int f = 0; f < cell->n_faces(); ++f)
        {
            auto neighbor = cell->neighbor(f);

            if(neighbor->at_boundary())

            continue;

            int neighbor_index = neighbor->active_cell_index();
            if (visited_cells[neighbor_index])
                                continue;
            // Check if any vertex of the current cell is close to any vertex of the neighbor cell
            for (unsigned int vertex1 = 0; vertex1 < cell->n_vertices(); ++vertex1)
            {
                for (unsigned int vertex2 = 0; vertex2 < neighbor->n_vertices(); ++vertex2)
                {
                    const dealii::Point<dim> &point1 = cell->vertex(vertex1);
                    const dealii::Point<dim> &point2 = neighbor->vertex(vertex2);

                    // Adjust the threshold based on your precision requirements
                    if (point1.distance(point2) < 1e-10)
                    {
                        int neighbor_index = neighbor->active_cell_index();
                        cell_contributions(cell_index, neighbor_index) += solid_3d.tangent_matrix(cell->vertex_dof_index(i, 0), neighbor->vertex_dof_index(i, 0));
                    }
                }
            }
        }
      }
  // Mark the cell as visited
    visited_cells[cell_index] = true;
  }
