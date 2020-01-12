run:
	python3 run.py
push:
	git push heroku master
all:
	git commit -a
upgit:
	git push -u origin master
rung:
	gunicorn run:server -b :8001 -n tcellbase
	#-w workers
	#-b :8001
	#-n dctoxo
	#--chdir  if we do not cd first...
	#--forwarded-allow-ips="130.239.192.39"
