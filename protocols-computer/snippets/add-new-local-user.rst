.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

Add a new local user
====================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Steps to add a new local user on our machines.


Steps
-----

1. Create new user in NFS interface by logging into machine via VPN. Fill out appropriate boxes (standard choices), add R/W to `homes`.  Determine the UID for this new user on the NFS by looking at ``/etc/passwd``.

2. On the local machine to which you are adding the user, run the following, makeing sure to set the UID correctly:

    .. code-block:: bash

        useradd -s /bin/zsh -g 100 -u <UID> jsalt
        passwd jsalt

3. Create or edit ``/usr/local/share/new_user.sh`` to contain (make sure you edit ``$NEWUSER``)

    .. code-block:: bash

        #!/bin/bash

        NEWUSER=jsalt
        mkdir -p /usr/local/ssh/users/$NEWUSER/.ssh
        touch /usr/local/ssh/users/$NEWUSER/.ssh/authorized_keys
        chown -R $NEWUSER:users /usr/local/ssh/users/$NEWUSER/
        chmod 0711 /usr/local/ssh/users/$NEWUSER/
        chmod 0700 /usr/local/ssh/users/$NEWUSER/.ssh
        chmod 0600 /usr/local/ssh/users/$NEWUSER/.ssh/authorized_keys

4. Make sure this is `chmod 0755`
5. Execute the file ``./new_user.sh``
6. Paste in the appropriate ``id_rsa.pub`` content.
