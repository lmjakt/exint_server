/* adds a suitable user log table that can be used to track usage */

/* the table has the following entries :
	1. user id (as specified on logging in)
	2. session_id specified as an int8. This should be a time specified by the connection object
	3. action. an int specifying the action (1 = login, 2 = logout ?)
	4. action time. The time when the action was taken.
	5. action date. The date when the action was taken.
	6. address. Text giving the ip address of the computer used

note that #2 should be the same for the login and logout action, but that we do not specify that these
have to be unique.
*/

drop table user_log;

begin;
create table user_log (
user_id int references users(index),
session_id int8,
action int,
time time,
date date,
address text,
socket int);
commit;

/* And that's all really */

