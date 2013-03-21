#ifndef CRYSTAL_TO_BOX_H
#define CRYSTAL_TO_BOX_H

#ifdef __cplusplus
extern "C" {
#endif

int lattice_to_cryst(double lattice[6], double matrix[9]);

#ifdef __cplusplus
}
#endif

#endif // CRYSTAL_TO_BOX_H
