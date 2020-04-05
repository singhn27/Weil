'''
Main Function: find_recurrence
Input: A list of data nums (say n_0, n_1, ..., n_l)
Output: The smallest list (a_0, a_1, ..., a_r) such that n_{k + r + 1} = a_{r}n_{k + r} + ... + a_0n_r 
	0 if no definite recurrence has been found

Algorithm Overview:
Take M_n,r to be the n x n matrix with (M_n,r)[i][j] = nums[i + j + r] (i.e. what make_matrix generates)
If there is a recurrence amongst the data of size n, then M_{n,r} will have determinant 0 for all r
The idea is we 
(1) Find the minimal n with that property - we don't check all r for time reasons, instead we check check_limit 
(2) Once we find n, we find an r such that M_{n-1,r} has nonzero determinant
(3) The desired value will be the final row of (M_{n-1,r+1})(M_{n-1, r})^{-1}
'''
from sage.all import *


data = [1, 324, 94538, 27342036, 7902286563, 2283769774008, 660009642608596, 1335199531221597380/7, 771745329973477544543/14, 1003654801709391177363266/63, 32228470855054446274072025/7, 922088779634264003152103052986/693, 177655771542871897226399473386043/462, 47675195263326471172200902896297624/429, 675128440123966217583531295758264620332/21021, 2926681787937393571125809628466744820624078/315315, 6766488293711253951946018417591857218404217201/2522520, 8310939246750847655828540501039686267396153273181/10720710, 43233505961597909605852099521954952446329307762692161/192972780, 13188621179729673423448808917121861288637394833136180999/203693490, 19602059250558217473051899100427301144562562239759316385203/1047566520, 11566031710298121559864372357840717084381875363789083732986423/2138781645, 73536829614075456710504696534588994912472160892551898904832727751/47053196190, 488799306444759560056180070980288189334125970777515647946396375195249/1082223512370, 3390311989500852299709405155619396384479683157934880741314170166437222959/25973364296880, 1224750206207182890857242946033205836040769584861966756734338938954667781521/32466705371100] 


def make_matrix(nums, start, dim):
	matrix_array = [[nums[i + j + start] for j in range(dim)] for i in range(dim)]
	return matrix(matrix_array)

def find_recurrence(nums, check_limit=4, start = 0, dim = 2):
	length = len(nums)
	final_dim = -1 #dimension of found recurrence
	
	#Part 1 in Algorithm Overview
	while (2*dim + start - 2 + check_limit < length and final_dim == -1): #While no out of bound errors AND no recurrence has been found
		M = make_matrix(nums, start, dim)		
		if (M.rank() != dim):
			final_dim = dim - 1
			for i in range(1,check_limit):
				M = make_matrix(nums, start + i, dim)
				if (M.rank() == dim):
					final_dim = -1
		dim += 1
	if (final_dim == -1):
		return 0

	#Part 2 in Algorithm Overivew
	M = make_matrix(nums, start, final_dim)
	N = make_matrix(nums, start + 1, final_dim)
	while (2 * final_dim + start  < length and M.rank() != final_dim): #While no out of bounds AND no nonzero det matrix has been found
		start += 1
		M = N
		N = make_matrix(nums, start + 1, final_dim)
	if (2*final_dim + start  >= length):
		 return 0

	#Part 3 in Algorithm Overview
	R = N*(~M)
	return R.row(-1)
		
if __name__ == '__main__':
	print(find_recurrence(data,2))
