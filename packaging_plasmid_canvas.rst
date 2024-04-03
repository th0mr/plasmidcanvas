Uploading a new version of plasmidcanvas to pypi
================================================

1 - Be a collaborator or owner of the plasmidcanvas pypi project (https://pypi.org/project/plasmidcanvas/). 
    - Contact thomrobinson76@gmail.com or penn.rainford@york.ac.uk to gain access.

2 - Setup poetry to use your pypi api key. See https://python-poetry.org/docs/repositories/#configuring-credentials

3 - Increment the version of the package in the pyproject.toml to the next major/minor/patch version.

3 - Run "poetry build" to build the package distribution.

4 - Run "poetry publish" to publish the package to pypi.

5 - Verify the new version is available at https://pypi.org/project/plasmidcanvas/.

Updating the pypi documentation in readthedocs
==============================================

TODO