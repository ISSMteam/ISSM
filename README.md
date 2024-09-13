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
Navigate to the parent directory where you want the ISSM repository to be located. If you plan to make contributions to the code base, we recommend that you check out a copy via SSH with,
```
git clone git@github.com:ISSMteam/ISSM.git
```
Note that you will you first need to <a href="https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account" target="_blank">add an SSH key to your GitHub account</a>.

If you plan only to use ISSM without making contributions, you also have the option of checking out a copy via HTTPS.

```
git clone https://github.com/ISSMteam/ISSM.git
```
Note that checkout via HTTPS does not require credentials, but does not allow commits without first setting up a personal access token.

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

If you find yourself making a lot of commits and pull requests, consider asking us to add you to the 'ISSM Contributors' group, which will allow you to make commits directly to the repository.

## Troubleshooting
### fatal: unable to access 'https[]()://github.com/ISSMteam/ISSM.git/': The requested URL returned error: 403
If you get this error on commit to the repository, it means you originally cloned via HTTPS. You have two options here,
1. <a href="https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens" target="_blank">Create a personal access token</a> and use it in place of your password when prompted for credentials.
2. Change the config for this clone of the repository to use the SSH protocol with,
````
git remote set-url origin git@github.com:ISSMteam/ISSM.git
````
Note that in the second case you will have to create an SSH key and add it to your GitHub account.

## Resources
 - Git Guides: https://github.com/git-guides
