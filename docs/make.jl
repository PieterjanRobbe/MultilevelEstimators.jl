using Documenter, MultilevelEstimators

makedocs(sitename = "MultilevelEstimators.jl",
		authors = "Pieterjam Robbe",
		pages = [
				 "Home" => "index.md",
				 "Manual" => Any[
								 "Guide" => "guide.md",
								 ],
				 "Reporter" => Any[
								   "Generating reports" => "reports.md",
								   ]
				 ]
		)
