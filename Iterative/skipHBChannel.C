bool skipHBChannel(int iphi,int ieta) {
  if ((iphi == -6) && (ieta<=-1) && (ieta>=-9))
    return true;
  else return false;
}
