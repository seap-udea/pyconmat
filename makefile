clean:
	find . -name "*~" -exec rm -rf {} \;

commit:
	git commit -am "Commit"
	git push 

pull:
	git pull
