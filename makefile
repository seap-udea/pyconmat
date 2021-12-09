clean:
	find . -name "*~" -exec rm -rf {} \;

commit:
	git commit -am "Commit"
	git push origin master

pull:
	git pull
