/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(jianopt@163.com)"
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include "pp_test.hpp"
#include "pp_mrc.hpp"

namespace deepem {

using namespace ::cppsci;
using namespace ::cppsci::mlearn;

struct Test {
#define TEST_PROPS \
    ELT(Str ,cnn_file         ,"KLH_rotation04.cnn"  ,"f"       ,"the name of the output cnn file" ) SEP \
    ELT(Str ,data_path        ,"KLHdata/"            ,"d"       ,"the path of the data set"        ) SEP \
    ELT(int ,boxsize          ,272                   ,"b"       ,"box size"                        ) SEP \
    ELT(int ,start_mic_num    ,57                    ,"start"   ,"start mic number"                ) SEP \
    ELT(int ,end_mic_num      ,57                    ,"end"     ,"end mic number"                  ) SEP \
    ELT(int ,dim_x            ,2048                  ,"x"       ,"dimension of x"                  ) SEP \
    ELT(int ,dim_y            ,2048                  ,"y"       ,"dimension of y"                  ) SEP \
    ELT(int ,scan_step        ,20                    ,"step"    ,"scan step"                       ) SEP \
    ELT(int ,range1           ,70                    ,"r1"      ,"range1"                          ) SEP \
    ELT(int ,range2           ,40                    ,"r2"      ,"range2"                          ) SEP \
    ELT(int ,min_std          ,22                    ,"min"     ,"min std"                         ) SEP \
    ELT(int ,max_std          ,34                    ,"max"     ,"max std"                         ) SEP \
    ELT(int ,name_length      ,2                     ,"nl"      ,"name length"                     ) SEP \
    ELT(Str ,name_prefix      ,""                    ,"np"      ,"name prefix"                     ) SEP \
    ELT(int ,rotation_angle   ,90                    ,"ra"      ,"rotation angle"                  ) SEP \
    ELT(int ,rotation_n       ,4                     ,"rn"      ,"rotation number"                 )
    CNN cnn;
    bool no_input = false;

#define ELT(type, name, value, abbr, intro) type name = value;
#define SEP
    TEST_PROPS
#undef ELT
#undef SEP

    Test(const Par &par) {
        read_par(par);
        STD_ ifstream ifile(cnn_file.c_str());
        if (!ifile) {
            std::cout << "\"" << cnn_file << "\" couldn't be found!" << std::endl;
            no_input = true;
        }
        else {
            ::cppsci::mlearn::parse(ifile, cnn);
            no_input = false;
        }
        ifile.close();
    }

    Test(const Par &par, CNN &cnn) {
        read_par(par);
        this->cnn = STD_ move(cnn);
        no_input = false;
    }

    void read_par(const Par &par) {
        MPI_LOG << "=========== Parameters ==============" << STD_ endl;
#define ELT(type, name, value, abbr, intro) par.set(name, #name); MPI_LOG << string::format("%-30s", #name) << name << std::endl;
#define SEP
        TEST_PROPS
#undef ELT
#undef SEP
        MPI_LOG << "=====================================" << STD_ endl;
    }

    Arrayd sub_img(const Arrayd &a, int row, int col, int l) {
        Arrayd p({l, l});
        int rows = a.shape(0);
        int cols = a.shape(1);
        const double *a_it = a.data()+row*cols+col;
        double *p_it = p.data();
        for (int i = 0; i < l; i++) {
            const double *a_it2 = a_it;
            for (int j = 0; j < l; j++) {
                *p_it = *a_it2;
                p_it++;
                a_it2++;
            }
            a_it += cols;
        }
        return STD_ move(p);
    }

    V<Vd> pick_box(const V<Vd> &box, double range) {
        V<Vd> newbox;
        for (int i = 0; i < box.size(); i++) {
            V<Vd> box_m;
            for (int j = 0; j < box.size(); j++) {
                if ((STD_ abs(box[j][0]-box[i][0])<=range) && (STD_ abs(box[j][1]-box[i][1])<=range)) {
                    box_m.push_back(box[j]);
                }
            }
            auto it = STD_ max_element(box_m.begin(), box_m.end(), [](const Vd &rt1, const Vd &rt2){
                return rt1[2]<rt2[2] || (rt1[2]==rt2[2]&&(rt1[0]<rt2[0]||rt1[1]<rt2[1]));
            });
            if (STD_ find_if(newbox.begin(), newbox.end(), [&it](const Vd &v){
                return eles_equal(*it, v);
            }) == newbox.end()) newbox.push_back(*it);
        }
        return STD_ move(newbox);
    }

    void run() {
        if (no_input) return;

        JN_INFOS("running");

        int size = boxsize*boxsize;
        JN_INFO(boxsize);

        for (auto && l : cnn.layers) l.a.clear();

#ifdef USEMPI
        // TODO : is this correct?
        // auto mpialloc = mpi_alloc(mpialloc.end-mpialloc.beg+1);
        auto mpialloc = mpi_alloc(end_mic_num-start_mic_num+1);
#pragma omp parallel for
        for (int image_num = start_mic_num+mpialloc.beg; image_num < start_mic_num+mpialloc.end; image_num++) {
#else
#pragma omp parallel for
        for (int image_num = start_mic_num; image_num <= end_mic_num; image_num++) {
#endif
#pragma omp critical
            MPI_LOG << "[*] Start to process Image " << image_num << std::endl;

            JN_INFO(image_num);
            Timer timer;

            Str c3 = string::format(string::merge("%0", name_length, 'd'), image_num);
            Str file = string::merge(data_path, "/mic/", name_prefix, c3, ".mrc");

            JN_INFO(c3);
            JN_INFO(file);

            if (!path_exists(file)) continue;

            auto && mig = mrc_read(file, 0);
            mig = image::rotate_nearest(mig,-90);
            JN_INFOA(mig);

            V<Vd> cal_r;
            JN_INFOS("calculating");
            JN_INFO(dim_x);
            JN_INFO(dim_y);
            JN_INFO(scan_step);
            //for (auto && l : cnn.layers) l.a.clear();
//#ifdef USEMPI
//            auto mpialloc = mpi_alloc(dim_x-boxsize+1, scan_step);
//            for (int i = mpialloc.beg; i < mpialloc.end; i += scan_step) {
//#else
//#pragma omp parallel for
            for (int i = 0; i < dim_x-boxsize+1; i += scan_step) {
//#endif
#pragma omp critical
                MPI_LOG << "[*] Image " << image_num << " column " << i+1 << std::endl;
                for (int j = 0; j < dim_y-boxsize+1; j += scan_step) {
                    JN_INFO(string::format("i:%d, j:%d", i, j));

                    auto &&par_img = sub_img(mig, i, j, boxsize);
                    JN_INFOA(par_img);

                    Arrayd cal({rotation_n, size});
                    double *cal_it = cal.data();

                    Arrayd par_std({rotation_n,1});
                    double *it = par_std.data();

                    for (int k = 0; k < rotation_n; k++) {
                        auto && im_rot = image::rotate_nearest(par_img, double(rotation_angle)*k);

                        double sum = 0;
                        for (int ll = 0; ll < size; ll++) sum += im_rot[ll];
                        double avg = sum/double(size);

                        sum = 0;
                        for (int ll = 0; ll < size; ll++) {
                            double x = im_rot[ll]-avg;
                            sum += x*x;

                            *cal_it = im_rot[ll];
                            cal_it++;
                        }

                        *it = STD_ sqrt(sum/double(size-1));
                        it++;
                    }

                    JN_INFO(par_std[0]);
                    JN_INFO(max_std);
                    JN_INFO(min_std);
                    if(par_std[0] < max_std && par_std[0] > min_std) {
                        cal = mapstd(cal);
                        cal.reshape({rotation_n,boxsize,boxsize});
                        JN_INFOA(cal);
                        CNN cnn2 = cnn;
                        auto && resu = cnn_test_m(cnn2,cal);
                        double sum = 0;
                        for (int k = 0; k < resu.size(); k++) sum += resu[k];
//#pragma omp critical
                        cal_r.push_back({double(i), double(j), sum/4.0, par_std[0]});
                    }
                }
            }

//#ifdef USEMPI
//            int size = cal_r.size();
//            int mpirank = mpi_rank();
//            int mpisize = mpi_size();
//            Vi sizes(mpisize);
//            MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
//            MPI_Bcast(sizes.data(), mpisize, MPI_INT, 0, MPI_COMM_WORLD);
//            double temp[4];
//            for (int i = 0; i < mpisize; i++) {
//                for (int j = 0; j < sizes[i]; j++) {
//                    if (i == mpirank) for (int k = 0; k < 4; k++) temp[k] = cal_r[j][k];
//                    MPI_Bcast(temp, 4, MPI_DOUBLE, mpirank, MPI_COMM_WORLD);
//                    if (i != mpirank) cal_r.push_back({temp[0], temp[1], temp[2], temp[3]});
//                }
//            }
//#endif

//#pragma omp critical
//            {
//                MPI_LOG << "[*] cal:" << std::endl;
//                for (auto && v : cal_r) {
//                    for (auto && i : v) {
//                        MPI_LOG << i << ' ';
//                    }
//                    MPI_LOG << std::endl;
//                }
//            }

            JN_INFOS("picking boxes");
            for (int ll = 0; ll < 9; ll++) {
                double threhold=0.1*(ll+1); //selection threhold
                JN_INFO(threhold);

                V<Vd> box;
                for (auto && p : cal_r) {
                    if (p[2]>threhold) box.push_back(p);
                }
                JN_INFO(box.size());
                if (box.size() < 2) continue;

                auto nnewbox = pick_box(pick_box(box, range1), range2);
                std::sort(nnewbox.begin(), nnewbox.end(), [](const Vd v1, const Vd v2){ return v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]); });

                JN_INFO(nnewbox.size());

                Str result_path = string::merge(data_path, string::format("/%02d", ll+1), "result");

                if (!path_exists(result_path)) path_make(result_path);

                Str box_name = string::merge(result_path, '/', name_prefix, c3, ".box");
                JN_INFO(box_name);

//#ifdef USEMPI
//                if (MPI_IS_ROOT) {
//#endif
                    STD_ ofstream fid(box_name.c_str());
                    for (int i = 0; i < nnewbox.size(); i++) {
                        Str line = string::format("%d %d %d %d", int(nnewbox[i][0])+1, int(nnewbox[i][1])+1, boxsize, boxsize);
                        JN_INFO(line);
                        fid << line << STD_ endl;
                    }
                    fid.close();
//#ifdef USEMPI
//                }
//#endif

            }
#pragma omp critical
            MPI_LOG << "[*] Elapsed time for processing Image " << image_num << ": " << timer.toc() << " seconds." << std::endl;
        }
    }
};

void pp_test_help() {
        MPI_LOG << "Partical picking procedure of ROME suite\n" << STD_ endl;
#define ELT(type, name, value, abbr, intro) MPI_LOG << string::format("-%-40s%s", string::merge(abbr,",-",JN_PP_STR2(name)).c_str(), intro) << STD_ endl;
#define SEP
    TEST_PROPS
#undef ELT
#undef SEP
        MPI_LOG << "\nAll the parameters could be put in a parameter file!" << STD_ endl;
}

void pp_test(mlearn::CNN &cnn, const Par &par) {
    Test(par, cnn).run();
}

void pp_test(const Par &par) {
    Test(par).run();
}

} // end namespace deepem

