## Resubmission
This is a resubmission. In this version I have:

* Improved the help pages. The most significant change was in the help page for the sdlrm() function,
  where I added the print() method and the See Also section. Additionally, I created a separate help 
  page for the summary.sdlrm() function.

* Added the main reference (Medeiros and Bourguignon, 2025) for the methods implemented in the package
  to all help pages. However, I have not included it in the DESCRIPTION file because there is no DOI
  or publicly available link yet. The paper will be submitted for publication this week, and I do 
  not have permission from my co-author to make a preprint available at this time.

* Included the \value{} section in the .Rd files for all exported methods and provided a detailed 
  explanation of the function outputs.
  
* Removed \dontrun{} from the examples, as they are all executable in under 5 seconds.

* Replaced print() and cat() with message() or warning() where appropriate. The only exceptions are 
  in the choose_mode() and envelope() functions, where logical arguments (trace and progressBar, 
  respectively) allow users to control text output and progress bars. These behaviors are documented.
  I understand that it is not necessary to make the change in the print() and summary() methods.
