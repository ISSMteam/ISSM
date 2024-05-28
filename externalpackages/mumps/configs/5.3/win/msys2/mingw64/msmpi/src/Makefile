#
#  This file is part of MUMPS 5.3.5, released
#  on Thu Oct 22 09:29:08 UTC 2020
#
#
topdir = ..
libdir = $(topdir)/lib
incdir = $(topdir)/include

default: d

.PHONY: default s d c z mumps_lib clean

s:
	$(MAKE) ARITH=s mumps_lib
d:
	$(MAKE) ARITH=d mumps_lib
c:
	$(MAKE) ARITH=c mumps_lib
z:
	$(MAKE) ARITH=z mumps_lib

include $(topdir)/Makefile.inc

mumps_lib:    $(incdir)/mumps_int_def.h \
              $(libdir)/libmumps_common$(PLAT)$(LIBEXT) \
              $(libdir)/lib$(ARITH)mumps$(PLAT)$(LIBEXT)

# Build $(incdir)/mumps_int_def.h, needed by mumps_c_types.h
$(incdir)/mumps_int_def.h: build_mumps_int_def 
	./build_mumps_int_def > $(incdir)/mumps_int_def.h
build_mumps_int_def:build_mumps_int_def.o
	$(CC) $(OPTL) build_mumps_int_def.o -o build_mumps_int_def

OBJS_COMMON_MOD = \
        ana_omp_m.o\
        ana_blk_m.o\
        ana_orderings_wrappers_m.o\
        double_linked_list.o\
        fac_asm_build_sort_index_ELT_m.o\
        fac_asm_build_sort_index_m.o\
        fac_descband_data_m.o\
        fac_future_niv2_mod.o\
        fac_ibct_data_m.o\
        fac_maprow_data_m.o\
        front_data_mgt_m.o\
        lr_common.o \
        mumps_comm_ibcast.o\
        mumps_l0_omp_m.o\
        mumps_memory_mod.o\
        mumps_mpitoomp_m.o\
        mumps_ooc_common.o\
        mumps_static_mapping.o\
        omp_tps_common_m.o


OBJS_COMMON_OTHER = \
        ana_blk.o\
        ana_orderings.o\
        ana_set_ordering.o\
        ana_AMDMF.o\
        bcast_errors.o\
        estim_flops.o\
        mumps_type_size.o \
        mumps_type2_blocking.o \
        mumps_version.o \
        mumps_print_defined.o \
        mumps_common.o\
        mumps_pord.o\
        mumps_metis.o\
        mumps_metis64.o\
        mumps_metis_int.o\
        mumps_scotch.o\
        mumps_scotch64.o\
        mumps_scotch_int.o\
        mumps_size.o\
        mumps_io.o\
        mumps_io_basic.o\
        mumps_io_thread.o\
        mumps_io_err.o\
        mumps_numa.o \
        mumps_thread.o \
        mumps_save_restore_C.o \
        mumps_config_file_C.o \
        mumps_thread_affinity.o \
        tools_common.o \
        sol_common.o

OBJS_MOD =   \
        $(ARITH)ana_aux.o \
        $(ARITH)ana_aux_par.o \
        $(ARITH)ana_lr.o\
        $(ARITH)fac_asm_master_ELT_m.o\
        $(ARITH)fac_asm_master_m.o\
        $(ARITH)fac_front_aux.o\
        $(ARITH)fac_front_LU_type1.o\
        $(ARITH)fac_front_LU_type2.o\
        $(ARITH)fac_front_LDLT_type1.o\
        $(ARITH)fac_front_LDLT_type2.o\
        $(ARITH)fac_front_type2_aux.o\
        $(ARITH)fac_sispointers_m.o\
        $(ARITH)fac_lr.o\
        $(ARITH)fac_mem_dynamic.o\
        $(ARITH)fac_omp_m.o\
        $(ARITH)fac_par_m.o\
        $(ARITH)fac_sol_l0omp_m.o\
        $(ARITH)lr_core.o\
        $(ARITH)lr_stats.o\
        $(ARITH)lr_type.o\
        $(ARITH)mumps_comm_buffer.o\
        $(ARITH)mumps_config_file.o\
        $(ARITH)mumps_load.o\
        $(ARITH)mumps_lr_data_m.o\
        $(ARITH)mumps_ooc_buffer.o\
        $(ARITH)mumps_ooc.o\
        $(ARITH)mumps_sol_es.o\
        $(ARITH)mumps_save_restore.o\
        $(ARITH)mumps_save_restore_files.o\
        $(ARITH)mumps_struc_def.o\
        $(ARITH)omp_tps_m.o\
        $(ARITH)sol_lr.o\
        $(ARITH)sol_omp_m.o\
        $(ARITH)static_ptr_m.o

OBJS_OTHER = \
        $(ARITH)ini_driver.o\
        $(ARITH)ana_driver.o\
        $(ARITH)fac_driver.o\
        $(ARITH)sol_driver.o\
        $(ARITH)sol_distrhs.o\
        $(ARITH)end_driver.o\
        $(ARITH)ana_aux_ELT.o\
        $(ARITH)ana_dist_m.o\
        $(ARITH)ana_LDLT_preprocess.o\
        $(ARITH)ana_reordertree.o\
        $(ARITH)arrowheads.o\
        $(ARITH)bcast_int.o\
        $(ARITH)fac_asm_ELT.o\
        $(ARITH)fac_asm.o\
        $(ARITH)fac_b.o\
        $(ARITH)fac_distrib_distentry.o\
        $(ARITH)fac_distrib_ELT.o\
        $(ARITH)fac_lastrtnelind.o\
        $(ARITH)fac_mem_alloc_cb.o\
        $(ARITH)fac_mem_compress_cb.o\
        $(ARITH)fac_mem_free_block_cb.o\
        $(ARITH)fac_mem_stack_aux.o\
        $(ARITH)fac_mem_stack.o\
        $(ARITH)fac_process_band.o\
        $(ARITH)fac_process_blfac_slave.o\
        $(ARITH)fac_process_blocfacto_LDLT.o\
        $(ARITH)fac_process_blocfacto.o\
        $(ARITH)fac_process_bf.o\
        $(ARITH)fac_process_end_facto_slave.o\
        $(ARITH)fac_process_contrib_type1.o\
        $(ARITH)fac_process_contrib_type2.o\
        $(ARITH)fac_process_contrib_type3.o\
        $(ARITH)fac_process_maprow.o\
        $(ARITH)fac_process_master2.o\
        $(ARITH)fac_process_message.o\
        $(ARITH)fac_process_root2slave.o\
        $(ARITH)fac_process_root2son.o\
        $(ARITH)fac_process_rtnelind.o\
        $(ARITH)fac_root_parallel.o\
        $(ARITH)fac_scalings.o\
        $(ARITH)fac_determinant.o\
        $(ARITH)fac_scalings_simScaleAbs.o\
        $(ARITH)fac_scalings_simScale_util.o\
        $(ARITH)fac_sol_pool.o\
        $(ARITH)fac_type3_symmetrize.o\
        $(ARITH)ini_defaults.o\
        $(ARITH)mumps_c.o\
        $(ARITH)mumps_driver.o\
        $(ARITH)mumps_f77.o\
        $(ARITH)mumps_gpu.o\
        $(ARITH)mumps_iXamax.o\
        $(ARITH)ana_mtrans.o\
        $(ARITH)ooc_panel_piv.o\
        $(ARITH)rank_revealing.o\
        $(ARITH)sol_aux.o\
        $(ARITH)sol_bwd_aux.o\
        $(ARITH)sol_bwd.o\
        $(ARITH)sol_c.o\
        $(ARITH)sol_fwd_aux.o\
        $(ARITH)sol_fwd.o\
        $(ARITH)sol_matvec.o\
        $(ARITH)sol_root_parallel.o\
        $(ARITH)tools.o\
        $(ARITH)type3_root.o

# NOTE:
# - We need to use Fortran compiler or otherwise link in Fortran libraries 
#	using C compiler (which does not make a lot of sense, but is possible)
# - In case of "Makefile:[line_number]: *** missing separator.  Stop.", ensure 
#	that $(FC) is preceded by a tab and not four spaces
#
$(libdir)/libmumps_common$(PLAT)$(LIBEXT):      $(OBJS_COMMON_MOD) $(OBJS_COMMON_OTHER)
	$(FC) -shared $(OUTF) $@ \
		$? $(LORDERINGS) $(LIBS) $(LIBLAPACK) $(LIBBLAS) $(LIBOTHERS) \
		-Wl,--out-implib=$(libdir)/libmumps_common$(PLAT).dll.a

$(libdir)/lib$(ARITH)mumps$(PLAT)$(LIBEXT):    $(OBJS_MOD) $(OBJS_OTHER)
	$(FC) -shared $(OUTF) $@ \
		$? -L$(libdir) -lmumps_common$(PLAT) $(LORDERINGS) $(LIBS) $(LIBLAPACK) $(LIBBLAS) $(LIBOTHERS) \
		-Wl,--out-implib=$(libdir)/lib$(ARITH)mumps$(PLAT).dll.a

# Dependencies between modules:
# i) arithmetic-dependent modules:
$(ARITH)ana_aux.o:              $(ARITH)mumps_struc_def.o \
                                mumps_static_mapping.o \
                                ana_orderings_wrappers_m.o \
                                ana_blk_m.o
$(ARITH)ana_aux_par.o:          $(ARITH)mumps_struc_def.o \
                                mumps_memory_mod.o \
                                ana_orderings_wrappers_m.o
$(ARITH)ana_lr.o:                $(ARITH)lr_core.o\
                                $(ARITH)lr_stats.o\
                                lr_common.o\
                                ana_orderings_wrappers_m.o \
                                ana_blk_m.o 
$(ARITH)fac_asm_master_ELT_m.o: omp_tps_common_m.o \
                                fac_ibct_data_m.o \
                                fac_asm_build_sort_index_ELT_m.o \
                                lr_common.o \
                                $(ARITH)fac_mem_dynamic.o \
                                $(ARITH)lr_core.o \
                                $(ARITH)ana_lr.o \
                                $(ARITH)mumps_lr_data_m.o \
                                $(ARITH)mumps_struc_def.o \
                                $(ARITH)omp_tps_m.o \
                                $(ARITH)mumps_comm_buffer.o \
                                $(ARITH)mumps_load.o
$(ARITH)fac_asm_master_m.o:     omp_tps_common_m.o \
                                fac_ibct_data_m.o \
                                fac_asm_build_sort_index_m.o \
                                lr_common.o \
                                $(ARITH)fac_mem_dynamic.o \
                                $(ARITH)lr_core.o \
                                $(ARITH)ana_lr.o \
                                $(ARITH)mumps_lr_data_m.o \
                                $(ARITH)mumps_struc_def.o \
                                $(ARITH)omp_tps_m.o \
                                $(ARITH)mumps_comm_buffer.o \
                                $(ARITH)mumps_load.o

$(ARITH)fac_front_aux.o:        $(ARITH)lr_type.o\
                                $(ARITH)lr_stats.o\
                                $(ARITH)mumps_comm_buffer.o\
                                $(ARITH)mumps_load.o\
                                $(ARITH)mumps_ooc.o\
                                mumps_ooc_common.o\
                                mumps_l0_omp_m.o

$(ARITH)fac_front_LU_type1.o   : $(ARITH)fac_front_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)fac_lr.o\
                                 $(ARITH)lr_type.o\
                                 $(ARITH)lr_stats.o\
                                 $(ARITH)ana_lr.o\
                                 $(ARITH)mumps_lr_data_m.o\
                                 mumps_l0_omp_m.o

$(ARITH)fac_front_LU_type2.o   : $(ARITH)fac_front_aux.o\
                                 $(ARITH)fac_front_type2_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)mumps_comm_buffer.o\
                                 mumps_comm_ibcast.o\
                                 $(ARITH)fac_lr.o\
                                 $(ARITH)lr_core.o\
                                 $(ARITH)lr_type.o\
                                 $(ARITH)lr_stats.o\
                                 $(ARITH)ana_lr.o\
                                 $(ARITH)mumps_lr_data_m.o\
                                 $(ARITH)mumps_struc_def.o

$(ARITH)fac_front_LDLT_type1.o : $(ARITH)fac_front_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)fac_lr.o\
                                 $(ARITH)lr_type.o\
                                 $(ARITH)lr_stats.o\
                                 $(ARITH)ana_lr.o\
                                 $(ARITH)mumps_lr_data_m.o\
                                 mumps_l0_omp_m.o

$(ARITH)fac_front_LDLT_type2.o : $(ARITH)fac_front_aux.o\
                                 $(ARITH)fac_front_type2_aux.o\
                                 $(ARITH)mumps_ooc.o\
                                 $(ARITH)mumps_comm_buffer.o\
                                 $(ARITH)mumps_load.o\
                                 $(ARITH)fac_lr.o\
                                 $(ARITH)lr_type.o\
                                 $(ARITH)lr_stats.o\
                                 $(ARITH)ana_lr.o\
                                 $(ARITH)mumps_lr_data_m.o\
                                 $(ARITH)mumps_struc_def.o

$(ARITH)fac_front_type2_aux.o : mumps_ooc_common.o\
                                $(ARITH)fac_front_aux.o\
                                $(ARITH)lr_type.o\
                                $(ARITH)mumps_struc_def.o\
                                $(ARITH)mumps_comm_buffer.o\
                                $(ARITH)mumps_load.o\
                                mumps_comm_ibcast.o\
                                fac_ibct_data_m.o

$(ARITH)fac_lr.o:               $(ARITH)lr_core.o\
                                $(ARITH)lr_type.o\
                                $(ARITH)mumps_lr_data_m.o\
                                $(ARITH)lr_stats.o

$(ARITH)fac_mem_dynamic.o:      $(ARITH)mumps_load.o\
                                $(ARITH)static_ptr_m.o

$(ARITH)fac_omp_m.o:            $(ARITH)fac_asm_master_m.o\
                                $(ARITH)fac_asm_master_ELT_m.o\
                                $(ARITH)fac_front_LU_type1.o\
                                $(ARITH)fac_front_LDLT_type1.o\
                                $(ARITH)mumps_load.o\
                                $(ARITH)omp_tps_m.o\
                                $(ARITH)lr_stats.o\
                                $(ARITH)mumps_struc_def.o\
                                omp_tps_common_m.o\
                                mumps_l0_omp_m.o

$(ARITH)fac_sol_l0omp_m.o:      $(ARITH)mumps_struc_def.o

$(ARITH)fac_omp_m.o:            $(ARITH)mumps_struc_def.o \
                                $(ARITH)fac_mem_dynamic.o\
                                $(ARITH)omp_tps_m.o\
                                omp_tps_common_m.o

$(ARITH)fac_par_m.o:            $(ARITH)mumps_load.o\
                                $(ARITH)mumps_ooc.o\
                                $(ARITH)fac_asm_master_m.o\
                                $(ARITH)fac_asm_master_ELT_m.o\
                                $(ARITH)omp_tps_m.o\
                                $(ARITH)fac_front_LU_type1.o\
                                $(ARITH)fac_front_LU_type2.o\
                                $(ARITH)fac_front_LDLT_type1.o\
                                $(ARITH)fac_front_LDLT_type2.o\
                                $(ARITH)fac_mem_dynamic.o\
                                $(ARITH)mumps_struc_def.o\
                                $(ARITH)lr_stats.o\
                                omp_tps_common_m.o\
                                mumps_l0_omp_m.o

$(ARITH)lr_core.o:              $(ARITH)lr_type.o\
                                $(ARITH)mumps_lr_data_m.o\
                                $(ARITH)lr_stats.o\
                                lr_common.o

$(ARITH)lr_stats.o:             $(ARITH)lr_type.o

$(ARITH)mumps_comm_buffer.o:    mumps_comm_ibcast.o \
                                $(ARITH)lr_type.o \
                                $(ARITH)lr_core.o \
                                $(ARITH)mumps_lr_data_m.o \
                                fac_ibct_data_m.o

$(ARITH)mumps_config_file.o:    $(ARITH)mumps_struc_def.o

$(ARITH)mumps_load.o:           $(ARITH)mumps_comm_buffer.o \
                                $(ARITH)mumps_struc_def.o \
                                fac_future_niv2_mod.o

$(ARITH)mumps_lr_data_m.o:      $(ARITH)lr_type.o\
                                front_data_mgt_m.o

$(ARITH)mumps_ooc_buffer.o:     mumps_ooc_common.o

$(ARITH)mumps_ooc.o:            $(ARITH)mumps_struc_def.o \
                                $(ARITH)mumps_ooc_buffer.o \
                                mumps_ooc_common.o

$(ARITH)mumps_sol_es.o:         $(ARITH)lr_type.o \
                                $(ARITH)mumps_lr_data_m.o

$(ARITH)mumps_save_restore.o:   $(ARITH)mumps_struc_def.o \
                                $(ARITH)mumps_save_restore_files.o \
                                $(ARITH)mumps_lr_data_m.o \
                                $(ARITH)mumps_ooc.o \
                                $(ARITH)fac_sol_l0omp_m.o \
                                front_data_mgt_m.o

$(ARITH)mumps_save_restore_files.o : $(ARITH)mumps_struc_def.o

$(ARITH)sol_lr.o:               $(ARITH)lr_type.o\
                                $(ARITH)lr_stats.o\
                                $(ARITH)mumps_lr_data_m.o

$(ARITH)sol_omp_m.o:            $(ARITH)mumps_struc_def.o





# Dependencies between modules:
# ii) arithmetic-independent modules:

ana_omp_m.o:                     double_linked_list.o
fac_asm_build_sort_index_ELT_m.o:omp_tps_common_m.o
fac_asm_build_sort_index_m.o:    omp_tps_common_m.o
fac_descband_data_m.o:           front_data_mgt_m.o
fac_ibct_data_m.o:               front_data_mgt_m.o
fac_maprow_data_m.o:             front_data_mgt_m.o
mumps_comm_ibcast.o:             fac_future_niv2_mod.o
mumps_static_mapping.o:          lr_common.o


# Compile modules before the rest
$(OBJS_COMMON_OTHER):$(OBJS_COMMON_MOD)
$(OBJS_OTHER):$(OBJS_COMMON_MOD) $(OBJS_MOD)

.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) -I. -I../include $(INCS) $(IORDERINGSF) $(ORDERINGSF) -c $*.F $(OUTF) $*.o
.c.o:
	$(CC) $(OPTC) -I../include $(INCS) $(CDEFS) $(IORDERINGSC) $(ORDERINGSC) -c $*.c $(OUTC) $*.o

$(ARITH)mumps_c.o:	mumps_c.c
	$(CC) $(OPTC) -I../include $(INCS) $(CDEFS) -DMUMPS_ARITH=MUMPS_ARITH_$(ARITH) \
	      $(IORDERINGSC) $(ORDERINGSC) -c mumps_c.c $(OUTC) $@

clean:
	$(RM) *.o *.mod $(incdir)/mumps_int_def.h
