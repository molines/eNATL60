# make file for tools

include make.macro

ALL: bat_mask_old_zone.exe

bat_mask_old_zone.exe: bat_mask_old_zone.f90
	$(F90) $(FFLAGS) bat_mask_old_zone.f90 -o  bat_mask_old_zone.exe

clean:
	\rm -f *~ *.exe
