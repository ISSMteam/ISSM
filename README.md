# Ice-sheet and Sea-level System Model - ISSM
[![C/C++ CI](https://github.com/ISSMteam/ISSM/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/ISSMteam/ISSM/actions/)

## Description
ISSM is a large-scale thermo-mechanical 2D/3D parallelized multi-purpose finite-element software dedicated to ice sheet and sea-level modeling.

## Contact
 - Forum:   https://issm.ess.uci.edu/forum/
 - GitHub:  https://github.com/ISSMteam/ISSM
 - Website: http://issm.jpl.nasa.gov/

## Checking Out a Copy of the Repository
Navigate to the parent directory where you want the ISSM repository to be located and run,
```
git clone https://github.com/ISSMteam/ISSM.git
```
or
```
git clone git@github.com:ISSMteam/ISSM.git
```

## Committing Changes to Repository
A good basic workflow for committing changes to the repository is,

1. Stash your local changes
```
git stash
```

2. Update your local branch
```
git pull
```

3. Merge your local changes
```
git stash apply
```

4. Add, commit, and push your changes
```
git add
git commit [-m <MSG>]
git push
```

You should then submit a pull request via GitHub so that project admins can review your changes and merge them into the main branch.
