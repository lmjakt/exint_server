create sequence protocol_type_seq;
create table protocol_types (
index int default nextval ('protocol_type_seq') unique,
type text,
user_id int references users (index),
entry_time timestamp default current_timestamp
);
create sequence protocol_step_seq;
create table protocol_step (
step int default nextval('protocol_step_seq') unique,
parent_step int default 1,
description text,
user_id int references users (index),
entry_time timestamp default current_timestamp
);
/* insert the root protocol step */
insert into protocol_step (description, user_id) values 
	('Root Inheritance Step', 1);
/* alter table so that steps must have some inheritance */
alter table protocol_step add constraint protocol_step_prnt_cons foreign key (parent_step) references protocol_step(step);
/* NOTE that for this to work the user interface will have to somehow encourage or force the user to indicate whether
   he/she is creating a new step or modifying an old one. I expect that this will be problematic... 
   but worth trying..
*/
create sequence protocol_seq;
create table protocols (
index int default nextval('protocol_seq') unique,
parent int default 1,
user_id int references users(index),
entry_time timestamp default current_timestamp,
type int references protocol_types(index),
name text not null, 
description text);
/* Insert a root protocol that all protocols will ultimately inherit, or be derived from */ 
insert into protocol_types (type, user_id) values (
	'Root Protocol Type', 1);
insert into protocols (user_id, type, name, description) values (
	1, 1, 'Root Protocol', 'Base inheritance protocol');
/* Now add a constraint into the protocol table such that every protocol
   has to inherit some protocol in the table, so that we can keep track of how 
   protocols are created by different users and so forth. 
*/
alter table protocols add constraint protocol_parent_const foreign key (parent) references protocols (index);
/* A protocol is made up of a sequeence of steps which are recorded in the protocol_step_series table
   these steps have fixed positions, -- rather than using a chain linked list type of arrangement. This
   is perhaps a little bit inflexible, but it is easier to implement and as any modification of a protocol
   will be interepreted as a creation of a new protocol with inheritance from the old protocol this 
   probably won't matter. The main trouble with this setup is that it does not track the specific change
   of the parent protocol, which would be preferrable, but more difficult to implement cleanly. */ 
create table protocol_step_series (
	protocol int references protocols (index),
	step int references protocol_step(step),
	step_number int);
/* I don't think there is an easy way to add a reasonable constraint on the step_number to ensure that
   there is always an appropriate series of events. One could make a linked list structure here, but 
   I fear we may run into issues of different types.. -- This is a bit of a kludge, but maybe it will
   do for now */


 
