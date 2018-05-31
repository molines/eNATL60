# make file for tools

include make.macro

ALL: bat_mask_old_zone.exe data_extend_zone.exe rnf_compute_runoff.exe rnf_reader.exe

bat_mask_old_zone.exe: bat_mask_old_zone.f90
	$(F90) $(FFLAGS) bat_mask_old_zone.f90 -o  bat_mask_old_zone.exe

data_extend_zone.exe: data_extend_zone.f90
	$(F90) $(FFLAGS) data_extend_zone.f90 -o  data_extend_zone.exe

rnf_reader.exe: rnf_reader.f90
	$(F90) $(FFLAGS) rnf_reader.f90 -o  rnf_reader.exe

rnf_compute_runoff.exe: rnf_compute_runoff.f90
	$(F90) $(FFLAGS) rnf_compute_runoff.f90 -o  rnf_compute_runoff.exe

clean:
	\rm -f *~ *.exe
