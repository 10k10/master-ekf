from matplotlib import pyplot as plt
import numpy as np
import ahrs

def Q2DCM(Q):
    # Extract the values from Q
    q0 = Q[0]
    q1 = Q[1]
    q2 = Q[2]
    q3 = Q[3]
     
    # First row of the rotation matrix
    r00 = 2 * (q0 * q0 + q1 * q1) - 1
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)
     
    # Second row of the rotation matrix
    r10 = 2 * (q1 * q2 + q0 * q3)
    r11 = 2 * (q0 * q0 + q2 * q2) - 1
    r12 = 2 * (q2 * q3 - q0 * q1)
     
    # Third row of the rotation matrix
    r20 = 2 * (q1 * q3 - q0 * q2)
    r21 = 2 * (q2 * q3 + q0 * q1)
    r22 = 2 * (q0 * q0 + q3 * q3) - 1
     
    # 3x3 rotation matrix
    rot_matrix = np.array([[r00, r01, r02],
                           [r10, r11, r12],
                           [r20, r21, r22]])
                            
    return rot_matrix

def Eul2DCM(theta1, theta2, theta3):
    r1 = np.array([[1, 0, 0],
                   [0, np.cos(theta1), np.sin(theta1)],
                   [0, -np.sin(theta1), np.cos(theta1)]])
    
    r2 = np.array([[np.cos(theta2), 0, -np.sin(theta2)],
                   [0, 1, 0],
                   [np.sin(theta2), 0, np.cos(theta2)]])

    r3 = np.array([[np.cos(theta3), np.sin(theta3), 0],
                   [-np.sin(theta3), np.cos(theta3), 0],
                   [0, 0, 1]])

    DCM = r3*r2*r1

    return DCM

lines = 0

time_vec = []

true_mag_x = []
true_mag_y = []
true_mag_z = []

true_sun_x = []
true_sun_y = []
true_sun_z = []

meas_mag_x = []
meas_mag_y = []
meas_mag_z = []

meas_ss_x1 = []
meas_ss_x1_valid = []
meas_ss_x2 = []
meas_ss_x2_valid = []
meas_ss_y1 = []
meas_ss_y1_valid = []
meas_ss_y2 = []
meas_ss_y2_valid = []
meas_ss_z1 = []
meas_ss_z1_valid = []
meas_ss_z2 = []
meas_ss_z2_valid = []

meas_sun_x = []
meas_sun_y = []
meas_sun_z = []

true_q1 = []
true_q2 = []
true_q3 = []
true_q4 = []

with open("./sim_data/data.42",'r') as data_file:
    for line in data_file:
        data = line.split(' ')

        if data[0] == 'TIME':
            time_vec.append(data[1])

        elif data[0] == 'SC[0].bvb':
            true_mag_x.append(float(data[2]))
            true_mag_y.append(float(data[3]))
            true_mag_z.append(float(data[4]))

        elif data[0] == 'SC[0].svb':
            true_sun_x.append(float(data[2]))
            true_sun_y.append(float(data[3]))
            true_sun_z.append(float(data[4]))

        elif data[0] == 'SC[0].AC.MAG[0].Field':
            meas_mag_x.append(float(data[2]))

        elif data[0] == 'SC[0].AC.MAG[1].Field':
            meas_mag_y.append(float(data[2]))

        elif data[0] == 'SC[0].AC.MAG[2].Field':
            meas_mag_z.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[0].Illum':
            meas_ss_x1.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[1].Illum':
            meas_ss_y1.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[2].Illum':
            meas_ss_z1.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[3].Illum':
            meas_ss_x2.append(-float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[4].Illum':
            meas_ss_y2.append(-float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[5].Illum':
            meas_ss_z2.append(-float(data[2]))
        
        elif data[0] == 'SC[0].AC.CSS[0].Valid':
            meas_ss_x1_valid.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[1].Valid':
            meas_ss_y1_valid.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[2].Valid':
            meas_ss_z1_valid.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[3].Valid':
            meas_ss_x2_valid.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[4].Valid':
            meas_ss_y2_valid.append(float(data[2]))

        elif data[0] == 'SC[0].AC.CSS[5].Valid':
            meas_ss_z2_valid.append(float(data[2]))
        
        elif data[0] == 'SC[0].B[0].qn':
            true_q1.append(float(data[2]))
            true_q2.append(float(data[3]))
            true_q3.append(float(data[4]))
            true_q4.append(float(data[5]))

        lines = lines + 1

for i in range(len(meas_ss_x1)):
    if meas_ss_x1_valid[i]:
        meas_sun_x.append(meas_ss_x1[i])
    elif meas_ss_x2_valid[i]:
        meas_sun_x.append(meas_ss_x2[i])
    else:
        meas_sun_x.append(0)

for i in range(len(meas_ss_y1)):
    if meas_ss_y1_valid[i]:
        meas_sun_y.append(meas_ss_y1[i])
    elif meas_ss_y2_valid[i]:
        meas_sun_y.append(meas_ss_y2[i])
    else:
        meas_sun_y.append(0)

for i in range(len(meas_ss_z1)):
    if meas_ss_z1_valid[i]:
        meas_sun_z.append(meas_ss_z1[i])
    elif meas_ss_z2_valid[i]:
        meas_sun_z.append(meas_ss_z2[i])
    else:
        meas_sun_z.append(0)

true_sun_x = np.reshape(np.transpose(np.array(true_sun_x)), (-1,1))
true_sun_y = np.reshape(np.transpose(np.array(true_sun_y)), (-1,1))
true_sun_z = np.reshape(np.transpose(np.array(true_sun_z)), (-1,1))
true_sun = np.concatenate((np.array(true_sun_x),np.array(true_sun_y),np.array(true_sun_z)),1)

meas_sun_x = np.reshape(np.transpose(np.array(meas_sun_x)), (-1,1))
meas_sun_y = np.reshape(np.transpose(np.array(meas_sun_y)), (-1,1))
meas_sun_z = np.reshape(np.transpose(np.array(meas_sun_z)), (-1,1))
meas_sun = np.concatenate((np.array(meas_sun_x),np.array(meas_sun_y),np.array(meas_sun_z)),1)

true_mag_x = np.reshape(np.transpose(np.array(true_mag_x)), (-1,1))
true_mag_y = np.reshape(np.transpose(np.array(true_mag_y)), (-1,1))
true_mag_z = np.reshape(np.transpose(np.array(true_mag_z)), (-1,1))
true_mag = np.concatenate((np.array(true_mag_x),np.array(true_mag_y),np.array(true_mag_z)),1)

meas_mag_x = np.reshape(np.transpose(np.array(meas_mag_x)), (-1,1))
meas_mag_y = np.reshape(np.transpose(np.array(meas_mag_y)), (-1,1))
meas_mag_z = np.reshape(np.transpose(np.array(meas_mag_z)), (-1,1))
meas_mag = np.concatenate((np.array(meas_mag_x),np.array(meas_mag_y),np.array(meas_mag_z)),1)

qn = np.zeros([len(true_mag_x),4])
qn = np.concatenate((np.reshape(true_q1,[-1,1]), np.reshape(true_q2,[-1,1]), np.reshape(true_q3,[-1,1]), np.reshape(true_q4,[-1,1])),1)

#sigma_noise = 0.05
#noise = np.random.normal(0,sigma_noise,[len(true_sun_x),3])

#meas_sun = meas_sun + noise

np.savetxt('./sim_data/true_sun.csv', true_sun, delimiter=',')
np.savetxt('./sim_data/meas_sun.csv', meas_sun, delimiter=',')
np.savetxt('./sim_data/true_mag.csv', true_mag, delimiter=',')
np.savetxt('./sim_data/meas_mag.csv', meas_mag, delimiter=',')
np.savetxt('./sim_data/qn.csv', qn, delimiter=',')

print('Done reading data\n')

B = np.zeros([3,3],float)
S = np.zeros([3,3],float)
K = np.zeros([4,4],float)
Z = np.zeros([3,1],float)
sigma = np.zeros([1,1],float)

q = ahrs.filters.davenport.Davenport(mag=meas_mag,acc=meas_sun)
print(q.Q.shape)

quat = np.zeros([len(time_vec),4],float)

quat[:,0] = q.Q[:,3]
quat[:,1] = q.Q[:,0]
quat[:,2] = q.Q[:,1]
quat[:,3] = q.Q[:,2]

""" for k in range(len(time_vec)):

    # isso esta errado, precisa voltar a rotaçao do quaternion para as medidas sem ruido
    b1 = meas_mag[k,:].reshape(3,1)
    b2 = meas_sun[k,:].reshape(3,1)
    n1 = true_mag[k,:].reshape(3,1)
    n2 = true_sun[k,:].reshape(3,1)

    

    B = b1*n1.T + b2*n2.T
    S = B + B.T
    sigma = np.trace(B)
    Z[0,0] = B[1,2]-B[2,1]
    Z[1,0] = B[2,0]-B[0,2]
    Z[2,0] = B[0,1]-B[1,0]

    K[0:3,0:3] = S-sigma*np.eye(3)
    K[0:3,3:4] = Z
    K[3:4,0:3] = Z.T
    K[3,3] = sigma

    eigs = np.linalg.eig(K)
    largest_ev = eigs[1][:,np.argmax(eigs[0])]
    #print(largest_ev)

    quat[k,:] = largest_ev """

print(type(quat))
fig, axs = plt.subplots(2,2)
axs[0,0].plot(quat[:,0])
axs[0,0].plot(qn[:,0])
axs[0,1].plot(quat[:,1])
axs[0,1].plot(qn[:,1])
axs[1,0].plot(quat[:,2])
axs[1,0].plot(qn[:,2])
axs[1,1].plot(quat[:,3])
axs[1,1].plot(qn[:,3])
#plt.show()