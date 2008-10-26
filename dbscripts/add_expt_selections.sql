/* create a set of sequences and tables that can hold selections of experiments
   and specify the order in which those should be displayed.
*/

--- create a sequence for the expt_selections. Data entry into the table must use this
--- but this does not have to use the default action..
begin;
create sequence expt_selection_seq;

--- make the table that will hold the expt_selections
create table expt_selections (
	selection int primary key default nextval('expt_selection_seq'),
	parent int references expt_selections(selection),
	selection_name text not null,
	selection_description text not null,
	user_id int references users(index),
	entry_time timestamp default current_timestamp,
	unique (selection, parent)
);

--- make the sequence for the actual lists...
create sequence expt_list_seq;

--- make the table that will hold the lists of experiments
create table expt_lists (
	list int primary key default nextval('expt_list_seq'),
	selection int,
	parent int,
	expt_id int references experiments (index),
	order_determinant float not null default 0,
	entry_time timestamp default current_timestamp,
	foreign key (selection, parent) references expt_selections (selection, parent),
	unique (parent, expt_id)
);

commit;

--- These tables need to be used with a great deal of care. In particular the syncrhonisation of the lists
--- will be very important as we could otherwise run into a great deal of trouble
