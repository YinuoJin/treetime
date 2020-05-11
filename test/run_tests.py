from test_treetime import *

import_short_test()

test_GTR()

test_ancestral()

test_seq_joint_reconstruction_correct()

test_seq_joint_reconstruction_ss_correct()

test_seq_joint_reconstruction_asvr_correct()

test_seq_joint_reconstruction_ss_asvr_correct()

test_seq_joint_reconstruction_mm_correct()

test_seq_joint_reconstruction_mm_asvr_correct()

# test is broken
#test_seq_joint_lh_is_max()

# test is broken
#test_seq_joint_lh_is_max_asvr()

print('\n\n TEST HAVE FINISHED SUCCESSFULLY\n\n')
