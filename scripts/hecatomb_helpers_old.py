def build_search_patterns(read_pattern_list, read_extension_list):
	
	search_pattern_list = []

	for curr_read_pattern in read_pattern_list:
		
		for curr_extension in read_extension_list:
			search_pattern = "*" + curr_read_pattern + "*" + curr_extension
			
			if not search_pattern in search_pattern_list:
				search_pattern_list.append(search_pattern)

	return(search_pattern_list)
