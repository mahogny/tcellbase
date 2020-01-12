run:
	python3 run.py
heroku:
	git push heroku master
all:
	git commit -a
github:
	git push -u origin master
rung:
	gunicorn run:server -b :8001 -n tcellbase
	#-w workers
	#-b :8001
	#-n dctoxo
	#--chdir  if we do not cd first...
	#--forwarded-allow-ips="130.239.192.39"

setheroku:
	heroku config:set ON_HEROKU=1
