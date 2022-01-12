#include "panhead.h"

void get_harm3d_data(double rr[],        double tt[],        double pp[],
		             double rho_ijk[],   double T_ijk[],     double bb_ijk[], 
		             double tau_ijk[],   double ut_ijk[],    double ur_ijk[], 
		             double uz_ijk[],    double up_ijk[],    long diskbody_ik[],
                     double sigtau_ik[], double Tdisk_ik[],  double emtop_ik[],
                     double embot_ik[],  double reftop_ik[], double refbot_ik[],
                     int rank) {
//  hid_t hfile_id, tfile_id, dataset_id;
    long i, j, k;

    /*
    hfile_id = H5Fopen(HFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    while (hfile_id < 0)
        hfile_id = H5Fopen(HFILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/r", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rr) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/th", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tt) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/phi", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pp) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/rho", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho_ijk) < 0)
        ;
    H5Dclose(dataset_id);
    */

    for (i = 0; i <= Nr; i++) {
        for (j = 0; j <= Nth; j++) {
            for (k = 0; k <= Nph; k++) {
                bb_ijk[indexijk(i,j,k)] = 0.0;
            }
        }
    }

    for (i = 0; i <= Nr; i++) {
        for (j = 0; j <= Nth+1; j++) {
            for (k = 0; k <= Nph; k++) {
                tau_ijk[indexthijk(i,j,k)] = 0.0;
            }
        }
    }

    /*
    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/ucon0", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ut_ijk) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/ucon1", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ur_ijk) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/ucon2", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uz_ijk) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/ucon3", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, up_ijk) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/diskbody_ik", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, diskbody_ik) < 0)
        ;
    H5Dclose(dataset_id);
    */

    for (i = 0; i <= Nr; i++) {
        for (k = 0; k <= Nph; k++) {
            sigtau_ik[indexr(i,k)] = 0.0;
        }
    }

    /*
    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/Tdisk_ik", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Tdisk_ik) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/emtop_ik", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, emtop_ik) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/embot_ik", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, embot_ik) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/emtop_ik", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, reftop_ik) < 0)
        ;
    H5Dclose(dataset_id);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(hfile_id, "/embot_ik", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, refbot_ik) < 0)
        ;
    H5Dclose(dataset_id);

    while (H5Fclose(hfile_id) < 0)
        ;
    */

    /*
    tfile_id = H5Fopen(TFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    while (tfile_id < 0)
        tfile_id = H5Fopen(HFILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    dataset_id = -1;
    while (dataset_id < 0)
        dataset_id = H5Dopen2(tfile_id, "/T_e", H5P_DEFAULT);
    while (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, T_ijk) < 0)
        ;
    H5Dclose(dataset_id);

    while (H5Fclose(tfile_id) < 0)
        ;
    */

    return;
}

