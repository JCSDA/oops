#!/usr/bin/env bash

# dirac_bump_cov
ln -sf parameters_bump_cov_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/dirac_bump_cov_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# dirac_bump_lct
ln -sf parameters_bump_lct_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/dirac_bump_lct_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# dirac_bump_hyb (hybrid localization is not very good, take localization alone instead)
ln -sf parameters_bump_loc_3d_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/dirac_bump_hyb_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# dirac_bump_loc_3d
ln -sf parameters_bump_loc_3d_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/dirac_bump_loc_3d_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# dirac_bump_loc_4d
ln -sf parameters_bump_loc_4d_nicas-2-sqrt_0001-0001_common.nc bump/dirac_bump_loc_4d_nicas-2-sqrt_0001-0001_common.nc

# 3densvar_bump
ln -sf parameters_bump_loc_3d_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/3densvar_bump_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# 3dvar_bump
ln -sf parameters_bump_cov_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/3dvar_bump_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# 3dvar_hybrid_bump (hybrid localization is not very good, take localization alone instead)
ln -sf parameters_bump_loc_3d_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/3dvar_hybrid_bump_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# 4densvar_bump
ln -sf parameters_bump_loc_4d_nicas-2-sqrt_0001-0001_common.nc bump/4densvar_bump_nicas-2-sqrt_0001-0001_common.nc

# 4densvar_advect_bump
ln -sf parameters_bump_loc_4d_nicas-2-sqrt_0001-0001_common.nc bump/4densvar_advect_bump_nicas-2-sqrt_0001-0001_common.nc

# 4dvar_drplanczos_bump
ln -sf parameters_bump_cov_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/4dvar_drplanczos_bump_nicas-2-sqrt_0001-0001_01_01_01_01.nc

# 4dvar_drplanczos_hybrid_bump (hybrid localization is not very good, take localization alone instead)
ln -sf parameters_bump_loc_3d_nicas-2-sqrt_0001-0001_01_01_01_01.nc bump/4dvar_drplanczos_hybrid_bump_nicas-2-sqrt_0001-0001_01_01_01_01.nc
