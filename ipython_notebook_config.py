# Ansible managed: /usr/lib/data-science-toolbox/bundles/base/ipython_notebook_config.py.j2 modified on 2014-03-13 13:47:46 by root on ip-10-181-106-120

c = get_config()
c.FileNotebookManager.notebook_dir = u'/home/ubuntu/notebooks'
c.IPKernelApp.pylab = 'inline'
c.NotebookApp.ip = '*'
c.NotebookApp.port = 8888
c.NotebookApp.open_browser = False
c.NotebookApp.password = u'sha1:b963e3827174:f033ba03720b12fd0e19c74571b20d6b142998c1'
c.NotebookApp.certfile = u'/home/ubuntu/.ssh/ipython_dst.pem'
c.NotebookManager.notebook_dir = u'/mnt/frontiers-review-2015'
