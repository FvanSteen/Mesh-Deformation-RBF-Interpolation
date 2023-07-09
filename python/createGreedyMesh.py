#%% creating the greedy mesh 2D

xElem = 200
yElem = 200

xNodes = xElem+1
yNodes = yElem+1

spacing = 0.04/8
blockL = int(0.2/spacing)
blockH = int(0.04/spacing)

nElem = int(xElem*yElem - (blockL*blockH))
nNodes = int(xNodes*yNodes - (blockL-1)*(blockH-1))

block_start_idx = int(((yElem/2)-blockH/2)*xNodes + (xElem-blockL)/2 )

o = open('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes\\400x400.su2', "w")
o.write("NDIME= 2\n")
o.write("NELEM= " + str(nElem) + "\n")

cnt = 0
skipped_idx = 0
for j in range( int(yElem/2-int(blockH)/2)):
    for i in range(xNodes-1):
        if (i < 81 or i > 120):
            o.write("9\t" + str(i + j*yNodes) + "\t" + str(i+1 + j*yNodes) + "\t" +  str(i+1+(1+j)*yNodes) + "\t" +  str(i+(1+j)*yNodes) + "\t" + str(cnt) +  "\n")
            cnt += 1            
        elif (j < 96 or j > 103):
            o.write("9\t" + str(i + j*yNodes) + "\t" + str(i+1 + j*yNodes) + "\t" +  str(i+1+(1+j)*yNodes) + "\t" +  str(i+(1+j)*yNodes) + "\t" + str(cnt) +  "\n")
            cnt += 1

idx = int(yElem/2-int(blockH)/2)*xNodes
for i in range(xNodes-1):
    if i < 80:
        o.write("9\t" + str(i + idx) + "\t" + str(i+1 + idx) + "\t" +  str(i+1+(1)*xNodes +idx) + "\t" +  str(i+(1)*xNodes +idx) + "\t" + str(cnt) +  "\n")
        cnt += 1
    elif i > 119:
        o.write("9\t" + str(i +idx) + "\t" + str(i+1 + idx) + "\t" +  str(i+1+(1)*xNodes + idx - (blockL-1)) + "\t" +  str(i+(1)*xNodes +idx - (blockL-1))  + "\t" + str(cnt) +  "\n")
        cnt += 1
idx = idx + xNodes
for j in range(blockH-2):
    for i in range(xNodes-1):
        if i < 80:
            o.write("9\t" + str(i +idx) + "\t" + str(i+1 + idx) + "\t" +  str(i+1+(1)*xNodes + idx - (blockL-1)) + "\t" +  str(i+(1)*xNodes +idx - (blockL-1))  + "\t" + str(cnt) +  "\n")
            cnt += 1
        elif i > 119:
            o.write("9\t" + str(i +idx- (blockL-1)) + "\t" + str(i+1 + idx- (blockL-1)) + "\t" +  str(i+1+(1)*xNodes + idx - 2*(blockL-1)) + "\t" +  str(i+(1)*xNodes +idx - 2*(blockL-1))  + "\t" + str(cnt) +  "\n")
            cnt += 1
    idx += (xNodes - (blockL-1))
    
for i in range(xNodes-1):
    if i < 80:
        o.write("9\t" + str(i + idx) + "\t" + str(i+1 + idx) + "\t" +  str(i+1+(1)*xNodes +idx- (blockL-1)) + "\t" +  str(i+(1)*xNodes +idx- (blockL-1)) + "\t" + str(cnt) +  "\n")
        cnt += 1
    elif i > 119:
        o.write("9\t" + str(i + idx - (blockL-1)) + "\t" + str(i+1 + idx - (blockL-1)) + "\t" +  str(i+1+(1)*xNodes +idx- (blockL-1)) + "\t" +  str(i+(1)*xNodes +idx- (blockL-1)) + "\t" + str(cnt) +  "\n")
        cnt += 1
idx += (xNodes - (blockL-1))

for j in range(int(yElem/2-int(blockH)/2)):
    for i in range(xNodes-1):
        o.write("9\t" + str(i + j*yNodes+idx) + "\t" + str(i+1 + j*yNodes+idx) + "\t" +  str(i+1+(1+j)*yNodes+idx) + "\t" +  str(i+(1+j)*yNodes+idx) + "\t" + str(cnt) +  "\n")
        cnt += 1

o.write("NPOIN= " + str(nNodes) + "\n")

cnt = 0
for j in range(yNodes):
    for i in range(xNodes):
        if (i < 81 or i > 119):
            o.write(str(i*spacing) + "\t" + str(j*spacing) + "\t" + str(cnt) + "\n")
            cnt += 1
        elif (j < 97 or j > 103):
            o.write(str(i*spacing) + "\t" + str(j*spacing) + "\t" + str(cnt) + "\n")
            cnt += 1

o.write("NMARK= 5\n")
o.write("MARKER_TAG= LOWER\n")
o.write("MARKER_ELEMS= "+str(xElem) + "\n")
for i in range(xElem):
    o.write("3\t" + str(i) + "\t" + str(i+1) + "\n")
    
o.write("MARKER_TAG= UPPER\n")
o.write("MARKER_ELEMS= "+str(xElem) + "\n")
s_idx = nNodes -xNodes
for i in range(xElem):
    o.write("3\t" + str(s_idx+i) + "\t" + str(s_idx+i+1) + "\n")
    
o.write("MARKER_TAG= LEFT\n")
o.write("MARKER_ELEMS= "+str(yElem) + "\n")
node = 0
for i in range(yElem):
    if(i < 97 or i > 103):
        o.write("3\t" + str(node) + "\t" + str(node+xNodes) + "\n")
        node += xNodes
    else:
        o.write("3\t" + str(node) + "\t" + str(node+xNodes-(int(blockL)-1)) + "\n")
        node += xNodes-(int(blockL)-1)
    
o.write("MARKER_TAG= RIGHT\n")
o.write("MARKER_ELEMS= "+str(yElem) + "\n")
node = xNodes-1
for i in range(yElem):
    if(i < 96 or i > 102):
        o.write("3\t" + str(node) + "\t" + str(node+xNodes) + "\n")
        node += xNodes
    else:
        o.write("3\t" + str(node) + "\t" + str(node+xNodes-(int(blockL)-1)) + "\n")
        node += xNodes-(int(blockL)-1)
      
o.write("MARKER_TAG= BLOCK\n")
o.write("MARKER_ELEMS= "+str(int(2*(blockL+blockH) )) + "\n")
for i in range(int(blockL)):
    o.write("3\t" + str(block_start_idx+i) + "\t" + str(block_start_idx+i+1) + "\n")
    
for i in range(int(blockH)):
    if i == int(blockH)-1:
        o.write("3\t" + str(block_start_idx+(int(blockL))+i*(xNodes-int(blockL)+1)) + "\t" + str(block_start_idx+(int(blockL))+i*(xNodes-int(blockL)+1)+xNodes) + "\n")
    else:
        o.write("3\t" + str(block_start_idx+(int(blockL))+i*(xNodes-int(blockL)+1)) + "\t" + str(block_start_idx+(int(blockL))+(i+1)*(xNodes-int(blockL)+1)) + "\n")
        
for i in range(int(blockL)):
    o.write("3\t" + str(block_start_idx+(int(blockL))+(int(blockH)-1)*(xNodes-int(blockL)+1)+xNodes - i) + "\t" + str(block_start_idx+(int(blockL))+(int(blockH)-1)*(xNodes-int(blockL)+1)+xNodes - i - 1)+ "\n")
    
for i in range(int(blockH)):
    if i == int(blockH)-1:    
        o.write("3\t" + str(block_start_idx+(int(blockL))+(int(blockH)-1)*(xNodes-int(blockL)+1)+xNodes - int(blockL) - (i)*(xNodes-int(blockL)+1)) + "\t" + str(block_start_idx) + "\n")
    else:
        o.write("3\t" + str(block_start_idx+(int(blockL))+(int(blockH)-1)*(xNodes-int(blockL)+1)+xNodes - int(blockL) - (i)*(xNodes-int(blockL)+1)) + "\t" + str(block_start_idx+(int(blockL))+(int(blockH)-1)*(xNodes-int(blockL)+1)+xNodes - int(blockL)  - (i+1)*(xNodes-int(blockL)+1)) + "\n")
o.close()