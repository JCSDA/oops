.PHONY: docs clean

docs:
	./makedocs
	@echo "Now open file://$(PWD)/html/index.html in your web browser"

clean:
	-rm *.pdf *.aux *.log
	-[ -d html ] && rm -r html

