# Ice-sheet and Sea-level System Model - ISSM
[![Ubuntu Basic](https://github.com/ISSMteam/ISSM/actions/workflows/ubuntu-basic.yml/badge.svg)](https://github.com/ISSMteam/ISSM/actions/workflows/ubuntu-basic.yml)
[![Ubuntu CodiPack](https://github.com/ISSMteam/ISSM/actions/workflows/ubuntu-codipack.yml/badge.svg)](https://github.com/ISSMteam/ISSM/actions/workflows/ubuntu-codipack.yml)

## Description
ISSM is a large-scale thermo-mechanical 2D/3D parallelized multi-purpose finite-element software dedicated to ice sheet and sea-level modeling.

## Contact
 - Forum:	https://issm.ess.uci.edu/forum
 - GitHub:	https://github.com/ISSMteam/ISSM
 - Website:	https://issm.jpl.nasa.gov

## Checking Out a Copy of the Repository
Navigate to the parent directory where you want the ISSM repository to be located and run,
```
git clone git@github.com:ISSMteam/ISSM.git
```

## Committing Changes to the Repository
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
git add [file]
git commit [-m "descriptive commit message"]
git push
```

5. If you have forked the ISSM repository, consider making sure that your commit passes CI workflows before submitting a pull request.

6. Submit a pull request via GitHub so that project admins can review your changes and merge them into the main branch.

## Troubleshooting
### fatal: unable to access 'https://github.com/ISSMteam/ISSM.git/': The requested URL returned error: 403
If you get this error on commit to the repository, it means you originally cloned via HTTPS. To fix this, switch to the SSH protocol with,
````
git remote set-url origin git@github.com:ISSMteam/ISSM.git
````

## Resources
 - Git Guides: https://github.com/git-guides
