/* make some tables and links for the sql queries.. 
   alter the ish-Image table so that it references the
   experiment tables.. and so on..
*/

create sequence ish_experiment_count;
create table ish_experiments (
	experiment int default nextval('ish_experiment_count') unique,
	user_id int references users (index),
	experiment_time timestamp,
	entry_time timestamp default current_timestamp,
	protocol int references protocols (index),
	free_comment text
);

/* insert a blank experiment into the table, -- for those data that we already have and 
   that we can't stricly define the experimental parameters for... */

insert into ish_experiments values (
	0, 1, '2003-07-04 12:50:46.591259+09', '2003-07-04 12:50:46.591259+09', 1,
	'This is a blank experiment for data that we can not map directly to a specific experiment. This is bad, but there is not so much we can do about it');

/* alter the image table to reference the experiment table.. */
alter table ish_images add column experiment int;
update ish_images set experiment = 0;
alter table ish_images add constraint experiment_cons foreign key (experiment) references ish_experiments (experiment);

/* alter the image table to reference the user who enters the table into the database.. */
alter table ish_images add column user_id int;
update ish_images set user_id=1;       /* the super user .. me ;-)  */
alter table ish_images add constraint user_cons foreign key (user_id) references users (index);

/* I think that may be all that we need to do at the moment.. just need user interfaces for these tables
   which respect the relationships.. */



