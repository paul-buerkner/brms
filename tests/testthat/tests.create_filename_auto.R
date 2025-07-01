test_that(" create_filename_auto works", {

  a  <- create_filename_auto( "a" , "b"  , F , nlist(a  = 5 ))
  a_exp <- nlist( file = "a" , file_refit = "b")
  expect_equal(a , a_exp )
  a  <- create_filename_auto( NA  , NULL   , F , nlist(a  = 5 , family = gaussian()  ))
  a_exp <- nlist( file = NA  , file_refit = NULL )
  expect_equal(a , a_exp )
  a  <- create_filename_auto( "a" , "b"  , F , nlist(a  = 5 , family = gaussian()  ))
  a_exp <- nlist( file = "a"  , file_refit = "b" )
  expect_equal(a , a_exp )

  a  <- create_filename_auto( "a" , "b"  , T , nlist(a  = 5 , family = gaussian()  ))
  a_exp <- nlist( file =  "cache-brm-result_a1d7ed8d4b70cd32.Rds"   , file_refit = "on_change" )
  expect_equal(a , a_exp )

  a  <- create_filename_auto( "a" , "b"  , T , nlist(a  = 5 , family = student()  ))
  a_exp <- nlist( file =  "cache-brm-result_a1d7ed8d4b70cd32.Rds"    , file_refit = "on_change" )
  expect_equal(a , a_exp )

  a  <- create_filename_auto( NULL , "something"  , T , nlist(a  = 5 , family = gaussian()  ))
  a_exp <- nlist( file =  "cache-brm-result_a1d7ed8d4b70cd32.Rds"  , file_refit = "on_change" )
  expect_equal(a , a_exp )
})
