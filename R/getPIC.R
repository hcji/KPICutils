f = 'E:/data/MetDIA_demo/data/30STD_mix 330ppb-1.mzML'
raw <- loadData(f)
mzs <- raw$mzs
ints <- raw$ints
scans <- raw$scans
times <- raw$times


getPIC <- function(mzs, scans, ints, times, noise, mztol=0.05, gap=1, min_pic_length=10){
  orders <- order(mzs)
  mzs <- mzs[orders]
  scans <- scans[orders]
  ints <- ints[orders]
  scantime <- times

  # set seeds
  seeds <- which(ints > noise)
  seeds <- seeds[order(-ints[seeds])]
  clu <- rep(0, length(ints))

  # detect pics
  clu <- getPIP(seeds,scans,mzs,clu,mztol,gap)
  orders <- order(clu)
  clu <- clu[orders]
  mzs <- mzs[orders]
  scans <- scans[orders]
  ints <- ints[orders]

  picind <- c(findInterval(0:max(clu),clu))
  pics <- lapply(1:(length(picind)-1),function(s){
    pic <- cbind(scans[(picind[s]+1):picind[s+1]],ints[(picind[s]+1):picind[s+1]],mzs[(picind[s]+1):picind[s+1]])
    pic <- pic[order(pic[,1]),]
    return(pic)
  })

  pic_lengths <- unlist(sapply(pics,length) / 3)
  pics <- pics[pic_lengths >= min_pic_length]
  return(pics)
}


getPeaks <- function(pics, noise, time_interval = 0.5, min_snr = 5, peakwidth=c(5,30)){
  pks <- lapply(seq_along(pics), function(i){
    pic <- pics[[i]]
    rtmin <- pic[1,1]
    rtmax <- pic[nrow(pic),1]
    rts <- seq(rtmin, rtmax, time_interval)
    intensities <- approx(pic[,1], pic[,2], rts)$y
    pks <- peak_detection(intensities, min_snr=min_snr, level=noise)
    keep <- pks$peakScale * time_interval > peakwidth[1] & pks$peakScale * time_interval < peakwidth[2]

    rt <- rts[pks$peakIndex[keep]]
    snr <- pks$snr[keep]
    intensity <- intensities[pks$peakIndex[keep]]
    mz <- rep(mean(pic[,3]), length(rt))
    profile <- rep(i, length(rt))

    cbind(rt, mz, intensity, snr, profile)
  })

  pks <- do.call(rbind, pks)
  return(pks)
}
