#include <vector>
#include <assert.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
#include <unistd.h>


inline int sub2ind(int i, int j, int nx) { return i * nx + j; }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
// void get_ordered_shoreline(int i_start, int j_start, std::vector<int> &ii,
//                            std::vector<int> &jj, bool *lake, int nx) {

  struct state {
  int i;
  int j;
  int dir;
  };
  int i_start;
  int j_start;
  std::vector<int> ii;
  std::vector<int> jj;
  bool *lake;
  int nx;

  // int sub2ind(const int row,const int col,const int cols,const int rows){
  //  return row*cols+col;
  //   }

  // Order for checking the points
  //       |      | 
  // -1, 1 | 0, 1 | 1, 1
  //       |      | 
  //---------------------
  //       |      | 
  // -1, 0 | 0, 0 | 1, 0
  //       |      | 
  //---------------------
  //       |      | 
  // -1,-1 | 0,-1 | 1,-1
  //       |      | 

  int dx[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  int dy[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  state start_state;
  start_state.i = i_start;
  start_state.j = j_start;
  start_state.dir = -1;
  // -1 indicates that the start_state dir hasn't been set yet. This will happen
  // in the first iteration of the loop.

  state next_state = start_state;
  next_state.dir = 0;

  while (true) {
    ii.push_back(next_state.i);
    jj.push_back(next_state.j);
    int dir =
        (next_state.dir + 5) % 8;  // start looking one after the last state
    int last_cell = 0;             // last cell is always land

    // Loop over all 8 directions.
    for (int i = 0; i < 8; i++) {
      int this_cell =
          lake[sub2ind(next_state.i + dx[dir], next_state.j + dy[dir], nx)];

      // Look for places where lake turns into land.
      if (last_cell && !this_cell) {
        // Advance the state.
        next_state.i += dx[dir];
        next_state.j += dy[dir];
        next_state.dir = dir;

        // Check the terminating condition.
        if ((next_state.i == start_state.i) &&
            (next_state.j == start_state.j) &&
            (next_state.dir == start_state.dir)) {
          return;
        }

        // The first time a state is found set the dir on start_state.
        if (start_state.dir == -1) {
          start_state.dir = dir;
        }

        // Print statement for debugging.
        printf("Adding state at position (%d,%d)\n", next_state.i, next_state.j);
        break;
      }

      // Prepare to look in the next direction.
      last_cell = this_cell;  // Move the cell information along.
      dir = (dir + 1) % 8;    // Increment the direction.
    }
  }

}

// int main() {
//   std::vector<int> ii;
//   std::vector<int> jj;

//   // Create a 7x7 example lake.
//   bool lake[49] = {0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 1, 1, 0, 0,
//                    0, 0, 1, 1, 1, 0, 0,
//                    0, 0, 1, 1, 1, 0, 0,
//                    0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 0, 0, 0, 0};
//   // Specify the lake width.
//   int nx = 7;

//   // Choose a starting coordinate.
//   int i_start = 2;
//   int j_start = 2;

//   // Get the ordered shoreline.
//   get_ordered_shoreline(i_start, j_start, ii, jj, lake, nx);

//   // Print the solution.
//   printf("shoreline indices:\n");
//   for (int i = 0; i < ii.size(); i++) {
//     printf("(%d, %d)\n", ii[i], jj[i]);
//   }

//   return 0;
// }