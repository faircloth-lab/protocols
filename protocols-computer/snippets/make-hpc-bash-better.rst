.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

Random Computer Snippets
========================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Our HPC systems use BASH (without the option of installing ZSH).  This is  a bit of a bummer, but you can also make BASH more ZSH-like with a few changes.

Steps
-----

1. Downloand and install `bash-it`_

.. _bash-it: https://github.com/Bash-it/bash-it


2. Because we use conda environments on the HPC, in `~/.bashrc` set
   
    .. code-block:: bash
   
        export BASH_IT_THEME='bobby-python'


3. While you are there, you may also want to add the following, which will give you prettier colors for `ls`:
   
   .. code-block:: bash
   
        eval "$(dircolors)"


4. Create `~/.inputrc` with the following contents - these changes let you use an anchor term and the up arrow to search backwards in history.  For example, if you type `cd` and hit the up arrow, you will search backwards in your history for all commands that start with `cd`.
   
    .. code-block:: bash
   
        "\e[A": history-search-backward
        "\e[B": history-search-forward

        set show-all-if-ambiguous on
        set completion-ignore-case on


