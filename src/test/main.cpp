#include "unit_test.h"
#include "cppsci_mpi.hpp"

int main(int argc, char **argv) {
    MPI_INITIALIZE(&argc, &argv);
    unit_test::run_test();
    MPI_FINALIZE;
}

