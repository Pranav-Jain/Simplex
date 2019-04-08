import numpy as np

basis_inds = []

def gen_basis(A,b,nv,ne):
	x_dash = np.zeros(ne)
	y = np.copy(b)

	x = np.concatenate((x_dash,y),axis=None)

	B = np.identity(nv)
	inds = np.arange(ne,ne+nv,1)

	c = np.zeros(ne)
	c = np.concatenate((c,np.ones(nv)),axis=None)

	A = np.concatenate((A,B),axis=1)

	simplex(A,c,b,x,A.shape[0],A.shape[1],B,inds)

	return(basis_inds)


def simplex(A,c,b,x,nv,ne,B,inds):

	B_inv = np.linalg.inv(B)

	x_b = np.dot(B_inv,b)

	for i in range(len(x_b)):
		if x_b[i]<0:
			print("INFEASIBLE")
			exit(0)

	k=0
	for i in inds:
		x[i] = x_b[k]
		k+=1

	
	arr = np.arange(0,ne,1)
	arr = np.delete(arr,inds,None)
	

	values = []
	
	for i in arr:
		c_b = []
		for j in inds:
			c_b.append(c[j])
		c_b = np.array(c_b)
		cj_ = c[i] - np.dot(c_b.T,np.dot(B_inv,A[:,i]))

		values.append(cj_)

	

	flag=0
	for i in range(len(values)):
		if values[i]<0:
			flag=1
	if flag==0:
		print('Optimal Values')
		print(x)
		print('Objective = ')
		print(np.dot(c.T,x))

		global basis_inds
		basis_inds = inds
		return(x)

	else:
		d = np.zeros(ne)
		index = np.argmin(values)
		
		d_b = np.dot(-1*B_inv,A[:,arr[index]])
		k=0
		for i in inds:
			d[i] = d_b[k]
			k+=1
		d[arr[index]]=1

		
		flag=0
		for i in d_b:
			if i<0:
				flag=1
		if flag==0:
			print('inf')
			return(-1)

		temp_arr = []
		index_arr=[]
		k=0
		for i in inds:
			if d[i]<0:
				temp_arr.append(-1*x[i]/d[i])
				index_arr.append(k)
			k+=1

		
		theta = np.min(np.array(temp_arr))

		index2 = np.argmin(np.array(temp_arr))
		index2 = index_arr[index2]

		x = x + theta*d

		B[:,index2] = A[:,arr[index]]

		inds = list(inds)
		inds[index2] = arr[index]

		# print('next iteration')

		simplex(A,np.copy(c),np.copy(b),np.copy(x),nv,ne,B,inds)
