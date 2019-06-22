# make file for tools

include make.macro

ALL: bat_mask_old_zone.exe data_extend_zone.exe rnf_compute_runoff.exe rnf_compute_runoff36.exe rnf_reader.exe rnf_mask.exe rnf_merge.exe

bat_mask_old_zone.exe: bat_mask_old_zone.f90
	$(F90) $(FFLAGS) bat_mask_old_zone.f90 -o  bat_mask_old_zone.exe

bat_mark_pool.exe: bat_mark_pool.f90
	$(F90) $(FFLAGS) bat_mark_pool.f90 -o  bat_mark_pool.exe

data_extend_zone.exe: data_extend_zone.f90
	$(F90) $(FFLAGS) data_extend_zone.f90 -o  data_extend_zone.exe

rnf_reader.exe: rnf_reader.f90
	$(F90) $(FFLAGS) rnf_reader.f90 -o  rnf_reader.exe

rnf_compute_runoff.exe: rnf_compute_runoff.f90
	$(F90) $(FFLAGS) rnf_compute_runoff.f90 -o  rnf_compute_runoff.exe

rnf_compute_runoff36.exe: rnf_compute_runoff36.f90
	$(F90) $(FFLAGS) rnf_compute_runoff36.f90 -o  rnf_compute_runoff36.exe

rnf_mask.exe: rnf_mask.f90
	$(F90) $(FFLAGS) rnf_mask.f90 -o  rnf_mask.exe

rnf_merge.exe: rnf_merge.f90
	$(F90) $(FFLAGS) rnf_merge.f90 -o  rnf_merge.exe

clean:
	\rm -f *~ *.exe
