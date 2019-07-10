.. include:: ../../links.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _SetupPersonalWebsite:

Setting Up Personal Websites
============================

:Author: Brant C. Faircloth
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

Modification History
--------------------

See `Setting Up Personal Websites`_.

.. _Setting Up Personal Websites: https://github.com/faircloth-lab/protocols/commits/master/protocols-data/lab/setting-up-personal-websites.rst


Purpose
-------

I've setup personal wordpress installations for everyone in the lab (that wants one) - primarily so that folks don't have to use weebly.


Steps
-----

#. Create appropriate entry in DNS
#. Log into webserver VPS
#. Log into MySQL as admin user
#. Create new database(s) for personal website(s), create a new user for that db, then assign privileges to DB user

    .. code-block:: sql

        CREATE DATABASE <name>_faircloth_lab;
        CREATE USER <name>_flab IDENTIFIED BY "STRONG PASSWORD";
        GRANT ALL PRIVILEGES ON <name>_faircloth_lab.* to <name>_flab@localhost IDENTIFIED BY "STRONG PASSWORD";

#. Copy over Nginx config for new user's website

    .. code-block:: bash

        cd /etc/nginx/conf.d/
        cp template.website.bak <name>.faircloth-lab.org

#. Edit file to reflect reasonable values (change ``<name>``)

    .. code-block:: text

        server {
            listen       80;
            server_name  <name>.faircloth-lab.org;
            access_log /var/www/html/<name>.faircloth-lab.org/logs/access.log;
            error_log /var/www/html/<name>.faircloth-lab.org/logs/error.log;

            # note that these lines are originally from the "location /" block
            root   /var/www/html/<name>.faircloth-lab.org/public_html;
            index index.php index.html index.htm;

            location / {
                try_files $uri $uri/ /index.php?$args;
            }
            error_page 404 /404.html;
            error_page 500 502 503 504 /50x.html;
            location = /50x.html {
                root /usr/share/nginx/html;
            }

            location ~ \.php$ {
                try_files $uri =404;
                fastcgi_pass unix:/var/run/php-fpm/php-fpm.sock;
                fastcgi_index index.php;
                fastcgi_param SCRIPT_FILENAME $document_root$fastcgi_script_name;
                include fastcgi_params;
            }
        }

#. Make appropriate directories in ``/var/www/html/``

    .. code-block:: bash

        cd /var/www/html
        mkdir -p <name>.faircloth-lab.org/{logs,public_html}

#. Download and unzip wordpress to ``public_html``
#. Change permissions of files within ``public_html``
#. Restart Nginx
#. Setup certbot for website and choose to redirect traffic from http to https:

    .. code-block:: bash
        
        certbot --nginx

#. Visit new site and setup

    .. code-block:: text

        https://<name>.faircloth-lab.org/

#. At new site, turn of ability to comment, remove sample comment, and turn off discussion for initial post
#. Copy over theme files
#. Adjust permissions

    .. code-block:: bash

        chown nginx:nginx <name>.faircloth-lab.org/{logs,public_html}
        cd public_html
        chown -R <admin>:<admin> -R *
        chown -R nginx:nginx wp-content
        find . -type d -exec chmod 755 {} \;
        find . -type f -exec chmod 644 {} \;