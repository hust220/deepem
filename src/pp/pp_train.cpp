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

#include "pp_train.hpp"
#include "pp_mrc.hpp"
#include "pp_box.hpp"
#include "../cppsci/cppsci_image.hpp"

namespace deepem {

using namespace ::cppsci;
using namespace ::cppsci::mlearn;

/**
 * Implementation of the training of partical picking.
 */
struct Train {
    bool debug_train = false;
#define TRAIN_PROPS \
    ELT(Str ,OUTPUT_CNN_NAME             ,"19s_rotation04"         ,"o"     ,"name of the output cnn file"           ) SEP \
    ELT(int ,boxsize                     ,160                      ,"bs"    ,"box size"                              ) SEP \
    \
    ELT(int ,dim_x                       ,1855                     ,"dx"    ,"x dimension"                           ) SEP \
    ELT(int ,dim_y                       ,1919                     ,"dy"    ,"y dimension"                           ) SEP \
    \
    ELT(int ,name_length                 ,5                        ,"nml"   ,"length of the name"                    ) SEP \
    ELT(Str ,name_prefix                 ,"image_"                 ,"nmp"   ,"prefix of the name"                    ) SEP \
    ELT(Str ,mic_path                    ,"19Sdata/mic/"           ,"mp"    ,"path of the mic"                       ) SEP \
    \
    ELT(int ,num_positive1               ,800                      ,"np1"   ,"number of positive1"                   ) SEP \
    ELT(int ,num_negative1               ,800                      ,"nn1"   ,"number of negative1"                   ) SEP \
    ELT(Str ,positive1_box_path          ,"19Sdata/positive/"      ,"p1bp"  ,"positive1 box path"                    ) SEP \
    ELT(Str ,negative1_box_path          ,"19Sdata/negative/"      ,"n1bp"  ,"negative1 box path"                    ) SEP \
    ELT(int ,positive1_mic_start_num     ,30001                    ,"p1msn" ,"positive1 mic start num"               ) SEP \
    ELT(int ,positive1_mic_end_num       ,30050                    ,"p1men" ,"positive1 mic end num"                 ) SEP \
    ELT(int ,negative1_mic_start_num     ,30001                    ,"n1msn" ,"negative1 mic start num"               ) SEP \
    ELT(int ,negative1_mic_end_num       ,30050                    ,"n1men" ,"negative1 mic end num"                 ) SEP \
    \
    ELT(int ,do_train_again              ,1                        ,"agn"   ,"whether to do train_agin (1:yes,0:no)" ) SEP \
    ELT(int ,num_positive2               ,800                      ,"np2"   ,"number of positive2"                   ) SEP \
    ELT(int ,num_negative2               ,800                      ,"nn2"   ,"number of negative2"                   ) SEP \
    ELT(Str ,positive2_box_path          ,"19Sdata/sel_positive/"  ,"p2bp"  ,"positive2 box path"                    ) SEP \
    ELT(Str ,negative2_box_path          ,"19Sdata/sel_negative/"  ,"n2bp"  ,"negative2 box path"                    ) SEP \
    ELT(int ,positive2_mic_start_num     ,30051                    ,"p2msn" ,"positive2 mic start num"               ) SEP \
    ELT(int ,positive2_mic_end_num       ,30090                    ,"p2men" ,"positive2 mic end num"                 ) SEP \
    ELT(int ,negative2_mic_start_num     ,30051                    ,"n2msn" ,"negative2 mic start num"               ) SEP \
    ELT(int ,negative2_mic_end_num       ,30100                    ,"n2men" ,"negative2 mic end num"                 ) SEP \
    \
    ELT(int ,rotation_angle              ,90                       ,"ra"    ,"rotation angle"                        ) SEP \
    ELT(int ,rotation_n                  ,360/rotation_angle       ,"rn"    ,"rotation number"                       ) SEP \
    \
    ELT(int ,num_p_test                  ,150                      ,"npt"   ,"number of positive test"               ) SEP \
    ELT(int ,num_n_test                  ,150                      ,"nnt"   ,"number of negative test"               ) SEP \
    \
    ELT(int ,FL_kernelsize               ,20                       ,"flk"   ,"FL kernel size"                        ) SEP \
    ELT(int ,TL_kernelsize               ,10                       ,"tlk"   ,"TL kernel size"                        ) SEP \
    ELT(int ,FIL_kernelsize              ,4                        ,"filk"  ,"FIL kernel size"                       ) SEP \
    \
    ELT(int ,SL_poolingsize              ,3                        ,"slp"   ,"SL pooling size"                       ) SEP \
    ELT(int ,FOL_poolingsize             ,2                        ,"folp"  ,"FOL pooling size"                      ) SEP \
    ELT(int ,SIL_poolingsize             ,2                        ,"silp"  ,"SIL pooling size"                      ) SEP \
    \
    ELT(int ,FL_feature_map              ,6                        ,"flfm"  ,"FL feature map"                        ) SEP \
    ELT(int ,TL_feature_map              ,12                       ,"tlfm"  ,"TL feature map"                        ) SEP \
    ELT(int ,FIL_feature_map             ,12                       ,"filfm" ,"FIL feature map"                       ) SEP \
    \
    ELT(int ,alpha                       ,1                        ,"a"     ,"alpha"                                 ) SEP \
    ELT(int ,batch_size                  ,50                       ,"bs"    ,"batch size"                            ) SEP \
    ELT(int ,num_epochs                  ,20                       ,"ne"    ,"number of epochs"						 )

#define ELT(type, name, value, abbr, intro) type name = value;
#define SEP
    TRAIN_PROPS
#undef ELT
#undef SEP

    Train(const Par &par) {
        MPI_LOG << "=========== Parameters ==============" << STD_ endl;
#define ELT(type, name, value, abbr, intro) par.set(name, #name, abbr); MPI_LOG << string::format("%-30s", #name) << name << std::endl;
#define SEP
            TRAIN_PROPS
#undef ELT
#undef SEP
        MPI_LOG << "=====================================\n" << STD_ endl;
    }

    void push_slis(L<Arrayd> &ls_train_x, const Arrayd &p, int beg, int slis, L<int> &ls_train_y, int y, int rotation_n) {
        int p_slis = p.shape(0);
        int rows = p.shape(1);
        int cols = p.shape(2);
        int size = rows*cols;

#ifdef RAND
        auto &&rands = randperm(p_slis);
#else
        auto &&rands = randperm(p_slis, p_slis*100);
#endif
        //for (auto && i : rands) std::cout << i << ' '; std::cout << std::endl;
        const double *p_it = p.data();
        for (int i = 0; i < slis; i++) {
            Arrayd temp({rows, cols});
            const double *p_it2 = p_it + rands[i+beg]*size;
            double *temp_it = temp.data();
            for (int j = 0; j < size; j++) {
                *temp_it = *p_it2;
                p_it2++;
                temp_it++;
            }
//            p_it += size;
            for (int j = 0; j < rotation_n; j++) {
                ls_train_x.push_back(image::rotate_nearest(temp,rotation_angle*j));
                ls_train_y.push_back(y);
            }
        }
    }

    /**
     * Set x array.
     */
    void set_x(const L<Arrayd> &ls, Arrayd &x, int rows, int cols) {
        int slis = ls.size();
        int size = rows*cols;
        x = Arrayd({slis, size});
        double *x_it = x.data();
        for (auto && l : ls) {
            const double *l_it = l.data();
            for (int i = 0; i < size; i++) {
                *x_it = *l_it;
                l_it++;
                x_it++;
            }
        }
        JN_INFOA(x);
        x = mapstd(x);
        x.reshape({slis, rows, cols});
    }

    /**
     * Set y array.
     */
    void set_y(const L<int> &ls, Arrayd &y) {
        int slis = ls.size();
        y = Arrayd({slis, 1});
        double *y_it = y.data();
        for (auto && l : ls) {
            *y_it = l;
            y_it++;
        }
    }

    /**
     * Create cnn.
     */
    CNN create_cnn() {
        Arrayd positive1, negative1, positive2, negative2;
        Arrayd ran1, ran2, ran3, ran4, a, train_x, train_y, test_x, test_y;
        CNN cnn;
        int i, j;
        double ri, p;
        Timer timer;

        MPI_LOG << "=========== Read data ==============" << STD_ endl;
        positive1=read_data(name_prefix, name_length, positive1_mic_start_num, positive1_mic_end_num, mic_path, positive1_box_path, boxsize,dim_x,dim_y);
        negative1=read_data(name_prefix, name_length, negative1_mic_start_num, negative1_mic_end_num, mic_path, negative1_box_path, boxsize,dim_x,dim_y);
        if (do_train_again) {
            positive2=read_data(name_prefix, name_length, positive2_mic_start_num, positive2_mic_end_num, mic_path, positive2_box_path, boxsize, dim_x, dim_y);
            negative2=read_data(name_prefix, name_length, negative2_mic_start_num, negative2_mic_end_num, mic_path, negative2_box_path, boxsize, dim_x, dim_y);
        }

        int rows = positive1.shape(1);
        int cols = positive1.shape(2);

        JN_INFO(do_train_again);
        JN_INFOA(positive1);
        JN_INFOA(negative1);
        if (do_train_again) {
            JN_INFOA(positive2);
            JN_INFOA(negative2);
        }
        MPI_LOG << "------------------------------------" << STD_ endl;
        MPI_LOG << "Elapsed time for reading data: " << timer.toc() << STD_ endl;
        MPI_LOG << "====================================\n" << STD_ endl;

        MPI_LOG << "=========== Prepare data ==============" << STD_ endl;
        L<Arrayd> ls_train_x;
        L<int> ls_train_y;

        push_slis(ls_train_x, positive1, 0, num_positive1, ls_train_y, 1, rotation_n);
        push_slis(ls_train_x, negative1, 0, num_negative1, ls_train_y, 0, rotation_n);

        if (do_train_again) {
            push_slis(ls_train_x, positive2, 0, num_positive2, ls_train_y, 1, rotation_n);
            push_slis(ls_train_x, negative2, 0, num_negative2, ls_train_y, 0, rotation_n);
        }

        set_x(ls_train_x, train_x, rows, cols);
        set_y(ls_train_y, train_y);

        JN_INFOA(train_x);
        JN_INFOA(train_y);

        L<Arrayd> ls_test_x;
        L<int> ls_test_y;

        push_slis(ls_test_x, positive1, num_positive1, num_p_test, ls_test_y, 1, 1);
        push_slis(ls_test_x, negative1, num_negative1, num_p_test, ls_test_y, 0, 1);

        set_x(ls_test_x, test_x, rows, cols);
        set_y(ls_test_y, test_y);

        JN_INFOA(test_x);
        JN_INFOA(test_y);

        MPI_LOG << "------------------------------------" << STD_ endl;
        MPI_LOG << "Elapsed time for preparing data: " << timer.toc() << STD_ endl;
        MPI_LOG << "====================================\n" << STD_ endl;

        MPI_LOG << "=========== CNN ==============" << STD_ endl;
        cnn.layers = {
            cnn_input_layer(),
            cnn_conv_layer(FL_feature_map, FL_kernelsize),
            cnn_samp_layer(SL_poolingsize),
            cnn_conv_layer(TL_feature_map, TL_kernelsize),
            cnn_samp_layer(FOL_poolingsize),
            cnn_conv_layer(FIL_feature_map, FIL_kernelsize),
            cnn_samp_layer(SIL_poolingsize)
        };

        CNNOpts opts;
        opts.alpha = alpha;
        opts.batchsize = batch_size;
        opts.numepochs = num_epochs;
        MPI_LOG << "[*] CNN init: " << timer.toc() << " seconds." << STD_ endl;

        cnn_setup(cnn, train_x, train_y);
        MPI_LOG << "[*] CNN setup: " << timer.toc() << " seconds." << STD_ endl;

        cnn_train(cnn, train_x, train_y, opts);
        MPI_LOG << "[*] CNN train: " << timer.toc() << " seconds." << STD_ endl;

#ifdef USEMPI
        if (MPI_IS_ROOT) {
#endif
            cnn_save(cnn, OUTPUT_CNN_NAME);
#ifdef USEMPI
        }
#endif
        MPI_LOG << "[*] Saving cnn file '" << OUTPUT_CNN_NAME << "': " << timer.toc() << " seconds." << STD_ endl;

        MPI_LOG << "------------------------------------" << STD_ endl;
        MPI_LOG << "Elapsed time for CNN: " << timer.toc() << STD_ endl;
        MPI_LOG << "====================================\n" << STD_ endl;

        return STD_ move(cnn);
    }

    /**
     * Read data.
     */
    Arrayd read_data(
        Str im, int b, int start_image_num, int end_image_num,
        Str mrc_directory, Str box_directory,
        int particle_image_size, int mic_width, int mic_length)
    {
        // im indicate the prefix  of the image name
        // b indicate the bites of the number in the image name
        int i;
        int num=0;
        int image_num;
        int boxsize=particle_image_size;

        L<Arrayd *> ls;

//        JN_INFO(start_image_num);
//        JN_INFO(end_image_num);
//        JN_INFO(mic_width);
//        JN_INFO(mic_length);

        for (image_num = start_image_num; image_num <= end_image_num; image_num++) {
            Str c = im;
            Str c3 = string::format(string::merge("%0", b, 'd'), image_num);
            Str file1 = string::merge(mrc_directory, c, c3, ".mrc");
            Str file2 = string::merge(box_directory, c, c3, ".box");
//            JN_INFO(c);
//            JN_INFO(c3);
            JN_INFO(file1);
            JN_INFO(file2);

            if (!path_exists(file1)) continue;
            if (!path_exists(file2)) continue;

            Arrayd mig=mrc_read(file1, 0);
//            JN_INFOA(mig);
            mig = image::rotate_nearest(mig,-90);
            mig = mapstd(mig);
//            JN_INFOA(mig);

            int rows = mig.shape(0);
            int cols = mig.shape(1);

            Arrayi box=box_read(file2);
//            JN_INFOA(box);
            if(box.shape(0) == 0) continue;

            int *box_it = box.data();
            for (i = 0; i < box.shape(0); i++) { //In case of some boxes out of bounds
                box_it[0]--;
                box_it[1]--;
                if(box_it[0]<0) box_it[0]=0; 
                if(box_it[1]<0) box_it[1]=0; 
                if(box_it[0]>mic_width-boxsize)  box_it[0]=mic_width-boxsize; 
                if(box_it[1]>mic_length-boxsize) box_it[1]=mic_length-boxsize; 

//                JN_INFO(box_it[0]);
//                JN_INFO(box_it[1]);
//                JN_INFO(mig(box_it[0], box_it[1]));

                // fetch the box from mig
                Arrayd *k = new Arrayd({boxsize, boxsize});
                double *k_it = k->data();
                double *mig_it = mig.data() + box_it[0]*cols+box_it[1];
                for (int row = 0; row < boxsize; row++) {
                    double *mig_it2 = mig_it;
                    for (int col = 0; col < boxsize; col++) {
                        *k_it = *mig_it2;
                        mig_it2++;
                        k_it++;
                    }
                    mig_it += cols;
                }
//                JN_INFOA(*k);
                ls.push_back(k);
                num++;
                box_it += 4;
            }
        }

        // move arrays in 'ls' to a and free spaces of these arrays at the mean time.
        Arrayd a({num,boxsize,boxsize});
        double *a_it = a.data();
        for (auto && k : ls) {
            double *k_it = k->data();
            for (int i = 0; i < k->size(); i++) {
                *a_it = *k_it;
                a_it++;
                k_it++;
            }
            delete k; // important !!!
        }

//        JN_INFOA(a);

        return STD_ move(a);
    }

};

void pp_train_help() {
    MPI_LOG << "Training procedure of the partical picking of ROME suite\n" << STD_ endl;
#define ELT(type, name, value, abbr, intro) MPI_LOG << string::format("-%-40s%s", string::merge(abbr,",-",JN_PP_STR2(name)).c_str(), intro) << STD_ endl;
#define SEP
    TRAIN_PROPS
#undef ELT
#undef SEP
        MPI_LOG << "\nAll the parameters could be put in a parameter file!" << STD_ endl;
}

cppsci::mlearn::CNN pp_train(const cppsci::Par &par) {
    auto && cnn = Train(par).create_cnn();
    return STD_ move(cnn);
}

} // end namespace deepem

