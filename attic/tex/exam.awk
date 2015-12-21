{
	if ($1 == "@include") {
		system("sed -f lf.sed "$2)
	} else {
		print $0
	}
}
