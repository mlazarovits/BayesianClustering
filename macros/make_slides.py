import os
import shutil
from pathlib import Path
from collections import defaultdict
import subprocess
import tempfile
import argparse
def escape_for_applescript(path):
    return f'"{str(path)}"'

def make_applescript_call_show(show):
    script = f'''
tell application "System Events"
    set visible of application process "Keynote" to {show}
end tell
'''
    with tempfile.NamedTemporaryFile("w", suffix=".applescript", delete=False) as f:
        f.write(script)
        script_path = f.name
    subprocess.run(["osascript", script_path])


def make_applescript_call_add_folder_title(folder_title):
    script = f'''
set slideTitle to "{folder_title}"
set imageFrames to {{ ¬
    {{{{80, 493}}, {{289, 237}}}}, ¬
    {{{{366, 493}}, {{289, 237}}}}, ¬
    {{{{655, 493}}, {{289, 237}}}}, ¬
    {{{{366, 257}}, {{289, 237}}}}, ¬
    {{{{80, 257}}, {{289, 237}}}}, ¬
    {{{{655, 257}}, {{289, 237}}}} ¬
}}

tell application "System Events"
	set visible of application process "Keynote" to false
end tell
tell application "Keynote"
    -- activate
    set theTheme to theme "ku-cms-minimalist" -- Change to your desired theme
    if (count of documents) = 0 then
        set thisDoc to make new document with properties {{document theme:theTheme}}
    else
        set thisDoc to front document
    end if

    tell thisDoc
        set thisSlide to make new slide with properties {{base slide:master slide "Plots"}}
        -- Iterate through text items and identify the first one (title box)
        repeat with ti in text items of thisSlide
            set tiPosition to position of ti
            -- Assuming the title box is the first text item
	    if (item 1 of tiPosition) is 75
                set object text of ti to slideTitle
                exit repeat
            end if
        end repeat
    end tell
end tell
'''

    with tempfile.NamedTemporaryFile("w", suffix=".applescript", delete=False) as f:
        f.write(script)
        script_path = f.name
    subprocess.run(["osascript", script_path])



def make_applescript_call_add_plots(pdf_paths, title_text):
	applescript_list = "{" + ", ".join([f'{escape_for_applescript(p)}' for p in pdf_paths]) + "}"
	#title_text = slide_title.replace("_", " vs ")

	#current slide list is a list of 1 string, needs to be a list of a list {{"pdf1","pdf2},{"pdf3"},{"pdf4"}}
	#list[i][j] = jth pdf on slide i
	script = f'''
	set pdfPaths to {applescript_list}
	set slideTitle to "{title_text}"
	set imageFrames_1 to {{ ¬
		{{{{83, 121}}, {{859, 582}}}} ¬
	}}
	set imageFrames_2 to {{ ¬
	    {{{{11, 118}}, {{510, 345}}}}, ¬
	    {{{{514, 408}}, {{510, 345}}}} ¬
	}}
	set imageFrames_3 to {{ ¬
	    {{{{262, 88}}, {{475,322}}}}, ¬
	    {{{{26, 428}}, {{475,322}}}}, ¬
	    {{{{538, 428}}, {{475,322}}}} ¬
	}}

	tell application "System Events"
    		set visible of application process "Keynote" to false
	end tell
	tell application "Keynote"
		-- activate
		set theTheme to theme "ku-cms-minimalist" -- Change to your desired theme
		if (count of documents) = 0 then
		    set thisDoc to make new document with properties {{document theme:theTheme}}
		else
		    set thisDoc to front document
		end if


		tell thisDoc
			set thisSlide to make new slide with properties {{base slide:master slide "Plots"}}
			-- Iterate through text items and identify the first one (title box)
			repeat with ti in text items of thisSlide
			    set tiPosition to position of ti
			    -- Assuming the title box is the first text item
			    if (item 1 of tiPosition) is 75
			        set object text of ti to slideTitle
			        exit repeat
			    end if
			end repeat
	
			-- will have different layouts based on # of plots in slide
			-- can add more #plot-dependent layouts here
			set pdfCount to count of pdfPaths
			if pdfCount is 1 then
				set frame to item 1 of imageFrames_1
				-- Get position and size from the frame
				set xpos to item 1 of item 1 of frame  -- X position
				set ypos to item 2 of item 1 of frame  -- Y position
				set imageWidth to item 1 of item 2 of frame  -- Width
				set imageHeight to item 2 of item 2 of frame  -- Height
			
				repeat with i from 1 to count of pdfPaths
					set thisPDF to item i of pdfPaths
					set pdfAlias to POSIX file thisPDF as alias
				
					tell thisSlide
					    -- Now apply these values to the new image
					    set newImg to make new image with properties {{file:pdfAlias}}
					    set position of newImg to {{xpos, ypos}}
					    set width of newImg to imageWidth
					    set height of newImg to imageHeight
					end tell
				end repeat
			else if pdfCount is 2 then
				repeat with i from 1 to count of pdfPaths
					set frame to item i of imageFrames_2
					-- Get position and size from the frame
					set xpos to item 1 of item 1 of frame  -- X position
					set ypos to item 2 of item 1 of frame  -- Y position
					set imageWidth to item 1 of item 2 of frame  -- Width
					set imageHeight to item 2 of item 2 of frame  -- Height
					set thisPDF to item i of pdfPaths
					set pdfAlias to POSIX file thisPDF as alias
				
					tell thisSlide
						-- Now apply these values to the new image
						set newImg to make new image with properties {{file:pdfAlias}}
						set position of newImg to {{xpos, ypos}}
						set width of newImg to imageWidth
						set height of newImg to imageHeight
					end tell
				end repeat
			else if pdfCount is 3 then
				repeat with i from 1 to count of pdfPaths
					set frame to item i of imageFrames_3
					-- Get position and size from the frame
					set xpos to item 1 of item 1 of frame  -- X position
					set ypos to item 2 of item 1 of frame  -- Y position
					set imageWidth to item 1 of item 2 of frame  -- Width
					set imageHeight to item 2 of item 2 of frame  -- Height
					set thisPDF to item i of pdfPaths
					set pdfAlias to POSIX file thisPDF as alias
				
					tell thisSlide
					    -- Now apply these values to the new image
					    set newImg to make new image with properties {{file:pdfAlias}}
					    set position of newImg to {{xpos, ypos}}
					    set width of newImg to imageWidth
					    set height of newImg to imageHeight
					end tell
				end repeat
			else
				return "Invalid # of pdfs = " & (pdfCount)
			end if
		end tell
	end tell
'''

	with tempfile.NamedTemporaryFile("w", suffix=".applescript", delete=False) as f:
	    f.write(script)
	    script_path = f.name

	# Print the AppleScript before running it for debugging
	# print("AppleScript:")
	# with open(script_path, 'r') as file:
	#     print(file.read())

	subprocess.run(["osascript", script_path])

# Function to sort suffix groups
def get_sorted_suffixes(suffix_groups):
    # Sort the suffixes by a predefined order of the prefixes
    sorted_suffixes = sorted(suffix_groups.items(), key=lambda x: x[0])  
    return sorted_suffixes


def clean_and_sort(groups):
	#remove empty slides
	newgroups = groups
	for i in groups:
		newgroups[i] = [x for x in newgroups[i] if len(x) > 0]
	newgroups = get_sorted_suffixes(newgroups)
	return newgroups

def make_slides(plot_dir_name):
	base_dir = os.getcwd()
	suffix_groups = defaultdict(list)
	# Define prefix (sample) order

	# Define plot directory and PNG files
	plot_dir = Path(plot_dir_name)
	full_dir = base_dir / plot_dir
	pdf_files = list(full_dir.glob("*.pdf"))
	pdf_files = sorted(pdf_files)
	#print("pdf_files",pdf_files)

	#maps slide title to plots on slide
	slideTitles_plots = {}
	slideTitles_plots["nJets"] = ["nJets"]
	slideTitles_plots["nSubclustersJet"] = ["nSubclustersJet"]
	slideTitles_plots["Jet Size"] = ["jetSize_ttbar"]
	slideTitles_plots["Jet Center"] = ["Jet_EtaCenter_ttbar","Jet_PhiCenter_ttbar","Jet_TimeCenter_ttbar"]
	slideTitles_plots["Subcluster Center"] = ["subClusterEtaCenter_ttbar","subClusterPhiCenter_ttbar","subClusterTimeCenter_ttbar"]
	slideTitles_plots["Jet Kinematics"] = ["Jet_mass","Jet_energy","Jet_pt"]
	slideTitles_plots["Gen Matching to Ws - dR, E ratio"] = ["genW_dR","genW_Eratio"]
	slideTitles_plots["Gen Matching to tops - dR, Eratio"] = ["genTop_dR","genTop_Eratio"]
	slideTitles_plots["Gen Matching to Ws - Kinematics"] = ["_genWMass","_genWE_","_genWPt"]
	slideTitles_plots["Gen Matching to tops - Kinematics"] = ["_genTopMass","_genTopE_","_genTopPt"]
	slideTitles_plots["Gen Matching to Ws - Center"] = ["_genWEtaCenter","_genWPhiCenter"]
	slideTitles_plots["Gen Matching to tops - Center"] = ["_genTopEtaCenter","_genTopPhiCenter"]
	slideTitles_plots["Gen Matching to Ws - subclusters"] = ["W_nSubclusters","W_subClusterEnergy"]

	for i in slideTitles_plots:
		suffix_groups[i].append([])

	# Process PDFs and group by suffix
	for pdf in pdf_files:
		stem = pdf.stem
		slidetag = ""
		slideidx = 0
		for i in slideTitles_plots:
			slide_check = [plot in stem for plot in slideTitles_plots[i]]
			if any(slide_check):
				#print("stem",stem,"title",i,"plots",slideTitles_plots[i],slide_check.index(True))
				slidetag = i
		#if no slidetag found
		if slidetag == "":
			print("Plot",stem,"does not have an associated slidetag yet.")
			continue
		slideidx = len(suffix_groups[slidetag]) - 1
		#if max # of plots per slide (based on slide tag) has been reached, add another slide to start to fill	
		if(len(suffix_groups[slidetag][slideidx]) >= len(slideTitles_plots[slidetag])):
			suffix_groups[slidetag].append([])
			slideidx += 1
		#do only ptstack hists for certain tags
		if "ptStack" not in stem:
			if "Center" in slidetag and "Gen Matching" not in slidetag:
				continue
			if slidetag == "Kinematics":
				continue
		if "Gen Matching" in slidetag and ("Kinematics" in slidetag or "Center" in slidetag):
			print(stem,slidetag)
			#only take matches of 'lead' aka high pt jets
			if "_lead" not in stem:
				continue
		#print(stem,"on slide",slidetag,"#",slideidx)
		suffix_groups[slidetag][slideidx].append(pdf)
		#for debugging
		#suffix_groups[slidetag][slideidx].append(stem)
	# Sorting the suffix groups by the predefined order
	#sorted_suffixes = get_sorted_suffixes(suffix_groups)
	sorted_suffixes = clean_and_sort(suffix_groups)
	
	#print("folder_title",plot_dir_name)	
	#make_applescript_call_add_folder_title(plot_dir_name)
	# Process plotting, iterating over sorted suffixes
	for suffix, slides in sorted_suffixes:
		#print("suffix",suffix,"slides",slides,"# slides",len(slides),"# pdfs on first slide",len(slides[0]))
		#print("suffix",suffix)#"slides",slides,"# slides",len(slides),"# pdfs on first slide",len(slides[0]))
		#if only 1 plot per slide
		#if len(pdfs) == 1:
		for pdflist in slides:
			make_applescript_call_add_plots(pdflist, suffix)
		#if 3 plots per slide (ie eta, phi, time centers or energy, mass pt)
		#break

def main():
	#plot_dir_names =[
	#"plots_local/condorSim_ttbar_defaultv9p8_defaultv7p2_bhcAlpha0p000_emAlpha0p000_thresh1p0_NperGeV0p250_beta0-2000_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NlnN/plots_5_30_25/"
	##"plots_local/testKeynoteMaker/"
	#]
	
	parser = argparse.ArgumentParser()
	parser.add_argument("--dirs","-d",help="dir(s) with plots to run over",required=True,nargs='+')
	args = parser.parse_args()
	make_applescript_call_show('false')
	for plot_dir_name in args.dirs:
		if(os.path.isdir(plot_dir_name)):
			print("making slides for",plot_dir_name)
			make_slides(plot_dir_name)
		else:
			print("dir",plot_dir_name,"doesn't exist")
		
	make_applescript_call_show('true')


if __name__ == "__main__":
	main()

