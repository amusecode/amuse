import numpy

class GridWithLevel(object):
    
    def __init__(self, grid, level, position,factor_per_level = 2):
        self.grid = grid
        self.level = level
        self.position = position
        self.factor_per_level = factor_per_level
        
    def position_at_level(self, level):
        return self.position_of_vector_at_level(self.position, level)
            
    def shape_at_level(self, level):
        return self.position_of_vector_at_level(self.grid.shape, level)
        
    def position_of_vector_at_level(self, vector, level):
        if level == self.level:
            return vector
        
        if level > self.level:
            dlevel = level - self.level
            factor = self.factor_per_level ** dlevel
            return tuple([x * factor for x in vector])
        else:
            dlevel = self.level - level 
            factor = self.factor_per_level ** dlevel
            return tuple([x / factor for x in vector])
            
    
    def position_at_level(self, level):
        return self.position_of_vector_at_level(self.position, level)
    
    def get_cell_for_position(self, position, level):
        box_of_one = tuple([1]*len(position))
        postion_on_my_level = self.position_of_vector_at_level(position, level)
        shape_on_my_level = self.position_of_vector_at_level(box_of_one, level)
        if shape_on_my_level == box_of_one:
            position_in_grid = tuple([x -y for x, y in zip(position, self.position)])
            return self.grid[position_in_grid]
        else:
            pass
            
class MultilevelBlockStructuredMesh(object):
    
    def __init__(self, shape, level = 1, position = (), grids_with_levels = [], factor_per_level = 2):
        self.level = level
        self._shape = shape
        self.factor_per_level = factor_per_level
        self.grids_with_levels = list(grids_with_levels)
        self.position = position
    
    @property
    def shape(self):
        return tuple(self._shape[len(self.position):])
        
    def add_grid(self, grid, level, *position):
        shape = grid.shape
        self.grids_with_levels.append(
            GridWithLevel(grid, level, position, self.factor_per_level)
        )
            
    def __getattr__(self, name):
        print name
    
    def __getitem__(self, index):
        print index
        if isinstance(index, int):
            newposition = list(self.position)
            newposition.append(index)
            if len(newposition) == len(self._shape):
                grids = self.select_grids_at_point(newposition)
                print len(grids)
                if len(grids) == 1:
                    grid = grids[0]
                    return grid.get_cell_for_position(newposition, self.level)
                elif len(grids) == 0:
                    raise Exception("Index out of range")
                else:
                    pass
            else:
                return MultilevelBlockStructuredMesh(
                    self._shape,
                    self.level,
                    newposition,
                    self.grids_with_levels,
                    self.factor_per_level
                )
    
    def to_grid(self):
        result = Grid(self.shape)
        grids_on_grid = self.select_grids_overlapping_rectangle(self.position)
        if len(grids_on_grid) == 1:
            return 
        
    def select_grids_at_point(self, indices):
        selected_grids = []
        for grid in self.grids_with_levels:
            position = grid.position_at_level(self.level)
            shape = grid.shape_at_level(self.level)
            
            must_select = True
            for index, pos, size in zip(indices, position, shape):
                if pos <= index and pos + size > index:
                    continue
                else:
                    must_select = False
                    break
            if must_select:
                selected_grids.append(grid)
        return selected_grids
            
    def select_grids_overlapping_rectangle(self, position, shape):
        selected_grids = []
        for grid in self.grids_with_levels:
            grid_position = grid.position_at_level(self.level)
            grid_shape = grid.shape_at_level(self.level)
            
            must_select = True
            for index, size, grid_index, grid_size in zip(position, shape, grid_position, grid_shape):
                if (
                    (grid_index           >= index and grid_index             < index+size )
                    or
                    (grid_index+grid_size >= index and grid_index + grid_size < index+size )
                   ):
                    continue
                else:
                    must_select = False
                    break
            if must_select:
                selected_grids.append(grid)
        return selected_grids
            
