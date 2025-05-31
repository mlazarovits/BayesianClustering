import os
import shutil
from pathlib import Path
from collections import defaultdict
import subprocess
import tempfile

def escape_for_applescript(path):
    return f'"{str(path)}"'

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

tell application "Keynote"
    activate
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
	print("title_text",title_text)
	#title_text = slide_title.replace("_", " vs ")


	print("pdf list",applescript_list)
	#current slide list is a list of 1 string, needs to be a list of a list {{"pdf1","pdf2},{"pdf3"},{"pdf4"}}
	#list[i][j] = jth pdf on slide i
	#TODO: modify structure (and revert applescript) to take in 1 slide (ie 1 pdf list) per add_plots call
	#NOT a list of slides
	#have catch in add_plots that checks how many pdfs are in the slide and formats them accordingly
	script = f'''
	set pdfPaths to {applescript_list}
	set slideTitle to "{title_text}"
	set imageFrames to {{ ¬
		{{{{83, 121}}, {{859, 582}}}} ¬
	}}

	tell application "Keynote"
		activate
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
			set pdfCount to count of pdfPaths
			-- return "# pdfs" & (pdfCount)
			if pdfCount is 1 then
				set frame to item 1 of imageFrames
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
			-- do formatting for 3 pdf slides
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

def make_slides(plot_dir_name):
	base_dir = os.getcwd()
	suffix_groups = defaultdict(list)
	# Define prefix (sample) order

	# Define plot directory and PNG files
	plot_dir = Path(plot_dir_name)
	full_dir = base_dir / plot_dir
	pdf_files = list(full_dir.glob("*.pdf"))
	#print("pdf_files",pdf_files)

	#maps slide title to plots on slide
	slideTitles_plots = {}
	slideTitles_plots["nJets"] = ["nJets"]
	slideTitles_plots["Center"] = ["EtaCenter_ttbar","PhiCenter_ttbar","TimeCenter_ttbar"]
	slideTitles_plots["Kinematics"] = ["Jet_mass","Jet_energy","Jet_pt"]

	for i in slideTitles_plots:
		suffix_groups[i].append([])

	print(suffix_groups)
	# Process PDFs and group by suffix
	for pdf in pdf_files:
		stem = pdf.stem
		slidetag = ""
		plotidx = 0
		slideidx = 0
		for i in slideTitles_plots:
			slide_check = [plot in stem for plot in slideTitles_plots[i]]
			if any(slide_check):
				#print("stem",stem,"title",i,"plots",slideTitles_plots[i],slide_check.index(True))
				slidetag = i
		#if max # of plots per slide (based on slide tag) has been reached, add another slide to start to fill	
		if(len(suffix_groups[slidetag][slideidx]) >= len(slideTitles_plots[slidetag])):
			suffix_groups[slidetag].append([])
			slideidx += 1
		print(stem,"on slide",slidetag,"slideidx",slideidx)
		suffix_groups[slidetag][slideidx].append(pdf)
		#if any of the single_plot tags are in stem, add pdf to single plot group with single_plot tag
		#if any of the triple_plot tags are in stem, add pdf to triple plot group
		#triple_check = [j in stem for i in triple_plot_slides for j in i]
		#triple_tag = [j for i in triple_plot_slides for j in i if j in stem]
		#single_tag = [j for i in single_plot_slides for j in i if j in stem]
		#if(any(single_tag)):
		#	tag = single_tag[0]
		#	print("single_tag",tag)
		#	suffix_groups[tag].append(pdf)
		#if(any(triple_check)):
		#	tag = triple_tag[0]
		#	print("triple_tag",tag)
		#	slide_idx = [idx for idx, i in enumerate(triple_plot_slides) if tag in i]
		#	print("slide_idx",slide_idx)
		#	suffix_groups[tag].append(pdf)
			
		#print("stem",stem,"pdf",pdf)
		#suffix_groups.append((stem,[pdf]))
		#suffix = stem

	print("suffix_groups",suffix_groups)
	# Sorting the suffix groups by the predefined order
	sorted_suffixes = get_sorted_suffixes(suffix_groups)

	print()
	print("sorted_suffixes",sorted_suffixes)
	print()
	#print("folder_title",plot_dir_name)	
	#make_applescript_call_add_folder_title(plot_dir_name)
	# Process plotting, iterating over sorted suffixes
	for suffix, slides in sorted_suffixes:
		#print("suffix",suffix,"slides",slides,"# slides",len(slides),"# pdfs on first slide",len(slides[0]))
		print("suffix",suffix)#"slides",slides,"# slides",len(slides),"# pdfs on first slide",len(slides[0]))
		#if only 1 plot per slide
		#if len(pdfs) == 1:
		for pdflist in slides:
			make_applescript_call_add_plots(pdflist, suffix)
		#if 3 plots per slide (ie eta, phi, time centers or energy, mass pt)
		#break

def main():
	plot_dir_names =[
	#"plots_local/condorSim_ttbar_defaultv9p8_defaultv7p2_bhcAlpha0p000_emAlpha0p000_thresh1p0_NperGeV0p250_beta0-2000_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NlnN/plots_5_30_25/"
	"plots_local/testKeynoteMaker/"
	]
	for plot_dir_name in plot_dir_names:
		if(os.path.isdir(plot_dir_name)):
			print("making slides for",plot_dir_name)
			make_slides(plot_dir_name)
		else:
			print("dir",plot_dir_name,"doesn't exist")
		


if __name__ == "__main__":
	main()

