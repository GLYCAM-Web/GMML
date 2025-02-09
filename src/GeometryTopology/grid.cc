#include "../../includes/GeometryTopology/coordinate.hpp"
#include "../../includes/GeometryTopology/cell.hpp"
#include "../../includes/GeometryTopology/grid.hpp"
#include "../../includes/common.hpp"

using GeometryTopology::Grid;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
Grid::Grid()
{
    min_corner_ = new GeometryTopology::Coordinate();
    max_corner_ = new GeometryTopology::Coordinate();
    cells_      = CellVector();
    assembly_   = NULL;
}

Grid::Grid(MolecularModeling::Assembly* assembly, GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max,
           double ion_radius, double ion_charge)
{
    min_corner_ = new GeometryTopology::Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new GeometryTopology::Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cells_      = CellVector();
    assembly_   = assembly;
    this->UpdateGrid(ion_charge);
}

Grid::Grid(MolecularModeling::Assembly* assembly, GeometryTopology::Coordinate* min, GeometryTopology::Coordinate* max,
           double cell_length, double cell_width, double cell_height)
{
    min_corner_ = new GeometryTopology::Coordinate(min->GetX(), min->GetY(), min->GetZ());
    max_corner_ = new GeometryTopology::Coordinate(max->GetX(), max->GetY(), max->GetZ());
    cells_      = CellVector();
    assembly_   = assembly;
    this->UpdateGrid(cell_length, cell_width, cell_height);
}

Grid::Grid(Grid& grid)
{
    //    std::cout<<"In Deep Copy"<<std::endl;
    GeometryTopology::Coordinate* tempMinCoordinate = new GeometryTopology::Coordinate(*grid.GetMinCorner());
    this->min_corner_                               = tempMinCoordinate;

    GeometryTopology::Coordinate* tempMaxCoordinate = new GeometryTopology::Coordinate(*grid.GetMaxCorner());
    this->max_corner_                               = tempMaxCoordinate;

    this->cells_ = grid.GetCells();

    MolecularModeling::Assembly* tempAssembly = new MolecularModeling::Assembly(*grid.GetAssembly());
    this->assembly_                           = tempAssembly;
}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
GeometryTopology::Coordinate* Grid::GetMinCorner()
{
    return min_corner_;
}

GeometryTopology::Coordinate* Grid::GetMaxCorner()
{
    return max_corner_;
}

Grid::CellVector Grid::GetCells()
{
    return cells_;
}

MolecularModeling::Assembly* Grid::GetAssembly()
{
    return assembly_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void Grid::SetMinCorner(GeometryTopology::Coordinate* min)
{
    min_corner_ = min;
}

void Grid::SetMaxCorner(GeometryTopology::Coordinate* max)
{
    max_corner_ = max;
}

void Grid::SetCells(CellVector cells)
{
    cells_.clear();
    for (CellVector::iterator it = cells.begin(); it != cells.end(); it++)
    {
        cells_.push_back(*it);
    }
}

void Grid::SetAssembly(MolecularModeling::Assembly* assembly)
{
    assembly_ = assembly;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void Grid::UpdateGrid(double length, double width, double height)
{
    double min_x = this->GetMinCorner()->GetX();
    double min_y = this->GetMinCorner()->GetY();
    double min_z = this->GetMinCorner()->GetZ();
    double max_x = this->GetMaxCorner()->GetX();
    double max_y = this->GetMaxCorner()->GetY();
    double max_z = this->GetMaxCorner()->GetZ();

    Grid::CellVector grid_cells = Grid::CellVector();
    for (double i = min_x; i <= max_x; i += length)
    {
        for (double j = min_y; j <= max_y; j += width)
        {
            for (double k = min_z; k <= max_z; k += height)
            {
                if (i + length <= max_x)
                {
                    if (j + width <= max_y)
                    {
                        if (k + height <= max_z)
                        {
                            Cell* new_cell =
                                new Cell(new GeometryTopology::Coordinate(i, j, k),
                                         new GeometryTopology::Coordinate(i + length, j + width, k + height));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                        else
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(i + length, j + width, max_z));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                    }
                    else
                    {
                        if (k + height <= max_z)
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(i + length, max_y, k + height));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                        else
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(i + length, max_y, max_z));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                    }
                }
                else
                {
                    if (j + width <= max_y)
                    {
                        if (k + height <= max_z)
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(max_x, j + width, k + height));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                        else
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(max_x, j + width, max_z));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                    }
                    else
                    {
                        if (k + height <= max_z)
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(max_x, max_y, k + height));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                        else
                        {
                            Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                      new GeometryTopology::Coordinate(max_x, max_y, max_z));
                            new_cell->SetGrid(this);
                            grid_cells.push_back(new_cell);
                        }
                    }
                }
            }
        }
    }
    this->SetCells(grid_cells);
    return;
}

void Grid::UpdateGrid(double ion_charge)
{
    double min_x = this->GetMinCorner()->GetX();
    double min_y = this->GetMinCorner()->GetY();
    double min_z = this->GetMinCorner()->GetZ();
    double max_x = this->GetMaxCorner()->GetX();
    double max_y = this->GetMaxCorner()->GetY();
    double max_z = this->GetMaxCorner()->GetZ();

    GeometryTopology::Coordinate* min_boundary = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* max_boundary = new GeometryTopology::Coordinate();
    this->assembly_->GetBoundary(min_boundary, max_boundary);

    if (max_x - min_x <= gmml::DEFAULT_BOX_LENGTH && max_y - min_y <= gmml::DEFAULT_BOX_WIDTH &&
        max_z - min_z <= gmml::DEFAULT_BOX_HEIGHT)
    {
        Grid::CellVector grid_cells = Grid::CellVector();
        for (double i = min_x; i <= max_x; i += gmml::DEFAULT_GRID_LENGTH)
        {
            for (double j = min_y; j <= max_y; j += gmml::DEFAULT_GRID_WIDTH)
            {
                for (double k = min_z; k <= max_z; k += gmml::DEFAULT_GRID_HEIGHT)
                {
                    if (i + gmml::DEFAULT_GRID_LENGTH <= max_x)
                    {
                        if (j + gmml::DEFAULT_GRID_WIDTH <= max_y)
                        {
                            if (k + gmml::DEFAULT_GRID_HEIGHT <= max_z)
                            {
                                Cell* new_cell =
                                    new Cell(new GeometryTopology::Coordinate(i, j, k),
                                             new GeometryTopology::Coordinate(i + gmml::DEFAULT_GRID_LENGTH,
                                                                              j + gmml::DEFAULT_GRID_WIDTH,
                                                                              k + gmml::DEFAULT_GRID_HEIGHT));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                            else
                            {
                                Cell* new_cell =
                                    new Cell(new GeometryTopology::Coordinate(i, j, k),
                                             new GeometryTopology::Coordinate(i + gmml::DEFAULT_GRID_LENGTH,
                                                                              j + gmml::DEFAULT_GRID_WIDTH, max_z));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                        }
                        else
                        {
                            if (k + gmml::DEFAULT_GRID_HEIGHT <= max_z)
                            {
                                Cell* new_cell =
                                    new Cell(new GeometryTopology::Coordinate(i, j, k),
                                             new GeometryTopology::Coordinate(i + gmml::DEFAULT_GRID_LENGTH, max_y,
                                                                              k + gmml::DEFAULT_GRID_HEIGHT));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                            else
                            {
                                Cell* new_cell = new Cell(
                                    new GeometryTopology::Coordinate(i, j, k),
                                    new GeometryTopology::Coordinate(i + gmml::DEFAULT_GRID_LENGTH, max_y, max_z));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                        }
                    }
                    else
                    {
                        if (j + gmml::DEFAULT_GRID_WIDTH <= max_y)
                        {
                            if (k + gmml::DEFAULT_GRID_HEIGHT <= max_z)
                            {
                                Cell* new_cell =
                                    new Cell(new GeometryTopology::Coordinate(i, j, k),
                                             new GeometryTopology::Coordinate(max_x, j + gmml::DEFAULT_GRID_WIDTH,
                                                                              k + gmml::DEFAULT_GRID_HEIGHT));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                            else
                            {
                                Cell* new_cell = new Cell(
                                    new GeometryTopology::Coordinate(i, j, k),
                                    new GeometryTopology::Coordinate(max_x, j + gmml::DEFAULT_GRID_WIDTH, max_z));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                        }
                        else
                        {
                            if (k + gmml::DEFAULT_GRID_HEIGHT <= max_z)
                            {
                                Cell* new_cell = new Cell(
                                    new GeometryTopology::Coordinate(i, j, k),
                                    new GeometryTopology::Coordinate(max_x, max_y, k + gmml::DEFAULT_GRID_HEIGHT));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                            else
                            {
                                Cell* new_cell = new Cell(new GeometryTopology::Coordinate(i, j, k),
                                                          new GeometryTopology::Coordinate(max_x, max_y, max_z));
                                if (!(new_cell->GetCellCenter()->GetX() > min_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetX() < max_boundary->GetX() &&
                                      new_cell->GetCellCenter()->GetY() > min_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetY() < max_boundary->GetY() &&
                                      new_cell->GetCellCenter()->GetZ() > min_boundary->GetZ() &&
                                      new_cell->GetCellCenter()->GetZ() < max_boundary->GetZ()))
                                {
                                    new_cell->SetGrid(this);
                                    grid_cells.push_back(new_cell);
                                }
                            }
                        }
                    }
                }
            }
        }
        this->SetCells(grid_cells);
        return;
    }
    else
    {
        Grid* temp_grid = new Grid();
        temp_grid->SetMaxCorner(this->GetMaxCorner());
        temp_grid->SetMinCorner(this->GetMinCorner());
        double min_x                = temp_grid->GetMinCorner()->GetX();
        double min_y                = temp_grid->GetMinCorner()->GetY();
        double min_z                = temp_grid->GetMinCorner()->GetZ();
        double max_x                = temp_grid->GetMaxCorner()->GetX();
        double max_y                = temp_grid->GetMaxCorner()->GetY();
        double max_z                = temp_grid->GetMaxCorner()->GetZ();
        Grid::CellVector grid_cells = Grid::CellVector();
        grid_cells.push_back(
            new Cell(temp_grid, temp_grid->GetMinCorner(),
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x) / 2,
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2)));
        grid_cells.push_back(new Cell(
            temp_grid,
            new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x) / 2,
                                             temp_grid->GetMinCorner()->GetY(), temp_grid->GetMinCorner()->GetZ()),
            new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x),
                                             temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                             temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2)));
        grid_cells.push_back(
            new Cell(temp_grid,
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX(),
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                                      temp_grid->GetMinCorner()->GetZ()),
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x) / 2,
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y),
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2)));
        grid_cells.push_back(
            new Cell(temp_grid,
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x) / 2,
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                                      temp_grid->GetMinCorner()->GetZ()),
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_y),
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y),
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2)));
        grid_cells.push_back(new Cell(
            temp_grid,
            new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX(), temp_grid->GetMinCorner()->GetY(),
                                             temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2),
            new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_y) / 2,
                                             temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                             temp_grid->GetMinCorner()->GetZ() + (max_z - min_z))));
        grid_cells.push_back(
            new Cell(temp_grid,
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x) / 2,
                                                      temp_grid->GetMinCorner()->GetY(),
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2),
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_y),
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z))));
        grid_cells.push_back(
            new Cell(temp_grid,
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX(),
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2),
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_y),
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y),
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z))));
        grid_cells.push_back(
            new Cell(temp_grid,
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_x) / 2,
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y) / 2,
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z) / 2),
                     new GeometryTopology::Coordinate(temp_grid->GetMinCorner()->GetX() + (max_x - min_y),
                                                      temp_grid->GetMinCorner()->GetY() + (max_y - min_y),
                                                      temp_grid->GetMinCorner()->GetZ() + (max_z - min_z))));

        temp_grid->SetCells(grid_cells);

        temp_grid->CalculateBoxCharge();
        temp_grid->CalculateBoxPotentialEnergy();
        Cell* best_box = temp_grid->GetBestBox(temp_grid, ion_charge);
        if (best_box != NULL)
        {
            Grid* new_grid = new Grid();
            new_grid->SetAssembly(this->assembly_);
            new_grid->SetMaxCorner(best_box->GetMaxCorner());
            new_grid->SetMinCorner(best_box->GetMinCorner());
            new_grid->UpdateGrid(ion_charge);
            this->SetCells(new_grid->GetCells());
            return;
        }
        else
        {
            //            std::cout << "Can't create an appropriate grid!" << std::endl;
            return;
        }
    }
}

void Grid::CalculateCellsCharge()
{
    for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
    {
        (*it)->CalculateCellCharge();
    }
}

void Grid::CalculateCellsPotentialEnergy(double ion_radius)
{
    for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
    {
        (*it)->CalculateCellPotentialEnergy(ion_radius);
    }
}

void Grid::CalculateBoxCharge()
{
    for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
    {
        (*it)->CalculateBoxCharge();
    }
}

void Grid::CalculateBoxPotentialEnergy()
{
    for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
    {
        (*it)->CalculateBoxPotentialEnergy();
    }
}

GeometryTopology::Cell* Grid::GetBestBox(Grid* grid, double ion_charge)
{
    if (ion_charge > 0)
    {
        double min_value = INFINITY;
        for (unsigned int i = 0; i < grid->GetCells().size(); i++)
        {
            if (grid->GetCells().at(i)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (min_value > grid->GetCells().at(i)->GetCellPotentialEnergy())
            {
                min_value = grid->GetCells().at(i)->GetCellPotentialEnergy();
            }
        }
        for (unsigned int i = 0; i < grid->GetCells().size(); i++)
        {
            if (grid->GetCells().at(i)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (fabs(min_value - grid->GetCells().at(i)->GetCellPotentialEnergy()) < gmml::THRESHOLD_PARTITIONING)
            {
                return grid->GetCells().at(i);
            }
        }
    }
    else if (ion_charge < 0)
    {
        double max_value = -INFINITY;
        for (unsigned int i = 0; i < grid->GetCells().size(); i++)
        {
            if (grid->GetCells().at(i)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (max_value < grid->GetCells().at(i)->GetCellPotentialEnergy())
            {
                max_value = grid->GetCells().at(i)->GetCellPotentialEnergy();
            }
        }
        for (unsigned int i = 0; i < grid->GetCells().size(); i++)
        {
            if (grid->GetCells().at(i)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (fabs(max_value - grid->GetCells().at(i)->GetCellPotentialEnergy()) < gmml::THRESHOLD_PARTITIONING)
            {
                return grid->GetCells().at(i);
            }
        }
    }
    return NULL;
}

std::vector<GeometryTopology::Coordinate*> Grid::GetBestPositions(double ion_charge)
{
    std::vector<GeometryTopology::Coordinate*> best_coordinates = std::vector<GeometryTopology::Coordinate*>();
    if (ion_charge > 0)
    {
        double min_value = INFINITY;
        for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
        {
            if ((*it)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (min_value > (*it)->GetCellPotentialEnergy())
            {
                min_value = (*it)->GetCellPotentialEnergy();
            }
        }
        for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
        {
            if ((*it)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (fabs(min_value - (*it)->GetCellPotentialEnergy()) < gmml::THRESHOLD)
            {
                best_coordinates.push_back(new GeometryTopology::Coordinate(
                    (*it)->GetCellCenter()->GetX(), (*it)->GetCellCenter()->GetY(), (*it)->GetCellCenter()->GetZ()));
            }
        }
    }
    else if (ion_charge < 0)
    {
        double max_value = -INFINITY;
        for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
        {
            if ((*it)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (max_value < (*it)->GetCellPotentialEnergy())
            {
                max_value = (*it)->GetCellPotentialEnergy();
            }
        }
        for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
        {
            if ((*it)->GetCellPotentialEnergy() == INFINITY)
            {
                continue;
            }
            if (fabs(max_value - (*it)->GetCellPotentialEnergy()) < gmml::THRESHOLD)
            {
                best_coordinates.push_back(new GeometryTopology::Coordinate(
                    (*it)->GetCellCenter()->GetX(), (*it)->GetCellCenter()->GetY(), (*it)->GetCellCenter()->GetZ()));
            }
        }
    }
    return best_coordinates;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void Grid::Print(std::ostream& out)
{
    for (CellVector::iterator it = cells_.begin(); it != cells_.end(); it++)
    {
        if ((*it)->GetCellPotentialEnergy() != INFINITY)
        {
            Cell* cell = *it;
            cell->GetCellCenter()->Print(out);
            out << ":  " << cell->GetCellPotentialEnergy() << std::endl;
        }
    }
}
