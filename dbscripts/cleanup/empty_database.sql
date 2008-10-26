/* sql code that cleans out the installation specific parts of an expression database
   primarily things like the expression data, user annotation, images and so forth
   this has to be done in a specific order and is a right pain to remembe by hand.. but.
*/

begin;    -- make sure we don't do something stupid.. 

-- first delete expression data and sample information
delete from data;
delete from cel_files;
delete from experiments;

-- in later versions we'll need to delete from expt_selections as well.. 

-- delete user annotation and comments. (remove the whole lot)
delete from user_annotation_genes;
delete from user_annotation;
delete from session_keywords;
delete from sessions where index != 0;

-- delete any images, image comments, and in_situ probes that might be lurking around..

-- the problem with ish_images is the braindead way in which postgresql handles images
delete from ish_annotation;
-- delete from ish_annotation_fields;   -- this is not strictly necessary, and it might be reasonable to include some of the fields.
-- I used to have a rule that removed any images associated with ish_images, but that seems
-- not to be there anymore
select lo_unlink(image) from ish_images;
select lo_unlink(thumbnail) from ish_images;
delete from ish_images;
delete from ish_tissue;

delete from ish_probe_text_annotation;
delete from ish_annotation_text_fields;
delete from ish_probe_classification;
delete from ish_probe_classes;
delete from ish_probe_num_annotation;
delete from ish_probe_num_fields;

-- lets see if we can delete from users and user privileges.
-- make sure to keep the database user with id = 1
-- make sure to set the password to the default affy_chip (-3946, 101151, -61370)
-- and set the user's name to chip and system administrator and stuff..
delete from user_privileges where user_index != 1;
delete from users where index != 1;
update users set user_name='chip';
update users set full_name='System Administrator';
update users set lab='Main Lab';
update users set key1 = -3946;
update users set key2 = 101151;
update users set key3 = -61370;

--- apart from the protocols and ish_experiments this is about all. But those tables are basically
--- empty anyway. 

-- seems to work so let's commit.
commit;
