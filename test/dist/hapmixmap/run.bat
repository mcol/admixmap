@echo ----------------------------------------------
@echo Step 1.
@echo    Running the initial training, from scratch.
@echo ----------------------------------------------
@pause
hapmixmap.exe training-initial.conf
@echo ----------------------------------------------
@echo Step 2.
@echo    Resuming the training, using saved state.
@echo ----------------------------------------------
@pause
hapmixmap.exe training-resume.conf
@echo ----------------------------------------------
@echo Step 3.
@echo    Testing and running the score test.
@echo ----------------------------------------------
@pause
hapmixmap.exe testing.conf
@echo ----------------------------------------------
@echo Step 4.
@echo    Post-processing the data.
@echo ----------------------------------------------
@pause
perl post-process.pl
