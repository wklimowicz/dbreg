# An unusable connection proves rejection precedes database setup/querying.
conn = structure(list(),class="dbreg_unusable_connection")
for (fun in list(dbreg::feols,dbreg::feglm)) {
  expect_error(do.call(fun,list(fml=y~i(f):x,conn=conn,table="absent")),"i.*terms must be additive")
  expect_error(do.call(fun,list(fml=y~i(f)*x,conn=conn,table="absent")),"i.*terms must be additive")
}
