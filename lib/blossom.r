blossom <- function(G, weighted = FALSE, maxcardinality = FALSE) {
  # A modified version without checking optimisation, otherwise for
  # weighted input the algorithm will fail
  mypath <- c()
  bestedgeto <- c()
  edges <- c(t(igraph::ends(G, igraph::E(G))))
  DEB <- FALSE
  DEBUG <- print
  CHECK_DELTA <- TRUE
  CHECK_OPTIMUM <- TRUE

  nedge <- length(edges) / 2
  # Sanity check
  if (weighted) {
    weights <- igraph::E(G)$weight
  } else {
    weights <- rep(1.0, nedge)
  }
  # TODO(bowendeng): stop if an edge is defined twice
  stopifnot(length(edges) %% 2 == 0 && nedge == length(weights))
  if (nedge == 0) {
    return(c())
  }

  # Count vertices.
  nvertex <- 0
  maxweight <- 0.0
  for (idx in 1:nedge) {
    i <- edges[2 * idx - 1]
    j <- edges[2 * idx]
    w <- weights[idx]
    stopifnot(i >= 0 && j >= 0 && i != j)
    nvertex <- max(nvertex, i, j)
    maxweight <- max(maxweight, w, na.rm = T)
  }

  # Construct a mapping from edge (v, w) where v < w to the corresponding weight
  edge_weight <- matrix(0, nrow = nvertex, ncol = nvertex)
  for (i in 1:nedge) {
    v <- edges[2 * i - 1]
    w <- edges[2 * i]
    if (v > w) {
      tmp <- v
      v <- w
      w <- tmp
    }
    wt <- weights[i]
    stopifnot(edge_weight[v, w] == 0 && edge_weight[w, v] == 0)
    edge_weight[v, w] <- wt
    edge_weight[w, v] <- wt
  }

  # If p is an edge endpoint,
  # endpoint[p] is the vertex to which endpoint p is attached.
  # Not modified by the algorithm.
  endpoint <- edges

  # If v is a vertex,
  # neighbend[v] is the list of remote endpoints of the edges attached to v.
  # Not modified by the algorithm.
  neighbend <- lapply(1:nvertex, function(x) list())
  for (idx in 1:nedge) {
    i <- edges[2 * idx - 1]
    j <- edges[2 * idx]
    neighbend[[i]][length(neighbend[[i]]) + 1] <- 2 * idx
    neighbend[[j]][length(neighbend[[j]]) + 1] <- 2 * idx - 1
  }

  # If v is a vertex,
  # mate[v] is the remote endpoint of its matched edge, or -1 of it is single
  # (i.e. endpoint[mate[v]] is v's partner vertex).
  # Initially all vertices are single; update during augmentation.
  mate <- rep(-1, nvertex)

  # If b is a top-level blossom,
  # label[b] is 0 if b is unlabeled (free);
  #             1 if b is an S-vertex/blossom;
  #             2 if b is a T-vertex/blossom.
  # The label of a vertex is found by looking at the label of its
  # top-level containing blossom.
  # If v is a vertex inside a T-blossom,
  # label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
  # Labels are assigned during a stage and reset after each augmentation.
  label <- rep(0, 2 * nvertex)

  # If b is a labeled top-level blossom,
  # labelend[b] is the remote endpoint of the edge through which b obtained
  # its label, or -1 if b's base vertex is single.
  # If v is a vertex inside a T-blossom and label[v] == 2,
  # labelend[v] is the remote endpoint of the edge through which v is
  # reachable from outside the blossom.
  labelend <- rep(-1, 2 * nvertex)

  # If v is a vertex,
  # inblossom[v] is the top-level blossom to which v belongs.
  # If v is a top-level vertex, v is itself a blossom (a trivial blossom)
  # and inblossom[v] == v.
  # Initially all vertices are top-level trivial blossoms.
  inblossom <- 1:nvertex

  # If b is a sub-blossom,
  # blossomparent[b] is its immediate parent (sub-)blossom.
  # If b is a top-level blossom, blossomparent[b] is -1.
  blossomparent <- rep(-1, 2 * nvertex)

  # If b is a non-trivial (sub-)blossom,
  # blossomchilds[b] is an ordered list of its sub-blossoms, starting with
  # the base and going round the blossom.
  blossomchilds <- lapply(1:(2 * nvertex), function(x) list())

  # If b is a (sub-)blossom,
  # blossombase[b] is its base VERTEX (i.e. recursive sub-blossom).
  blossombase <- c(1:nvertex, rep(-1, nvertex))

  # If b is a non-trivial (sub-)blossom,
  # blossomendps[[b]] is a list of endpoints on its connecting edges,
  # such that blossomendps[[b]][i] is the local endpoint of blossomchilds[[b]][i]
  # on the edge that connects it to blossomchilds[b][wrap(i+1)].
  blossomendps <- lapply(1:(2 * nvertex), function(x) list())

  # If v is a free vertex (or an unreachable vertex inside a T-blossom),
  # bestedge[v] is the edge to an S-vertex with lease slack,
  # or -1 if there is no such edge.
  # If b is a (possibly trivial) top-level S-blossom,
  # bestedge[b] is the lease-slack edge to a different S-blossom,
  # or -1 if there is no such edge.
  # This is used for efficient computation of delta2 and delta3.
  bestedge <- rep(-1, 2 * nvertex)

  # If b is a non-trivial top-level S-blossom,
  # blossombestedges[[b]] is a list of lease-slack edges to neighbouring
  # S-blossoms, or None if no such list has been computed yet.
  # This is used for efficient computation of delta3.
  blossombestedges <- lapply(1:(2 * nvertex), function(x) list())

  # List of currently unused blossom numbers.
  unusedblossoms <- (nvertex + 1):(2 * nvertex)

  # If v is a vertex,
  # dualvar[v] = 2 * u(v) where u(v) is the v's variable in the dual
  # optimization problem (multiplication by two ensures integer values
  # throughout the algorithm if all edge weights are integers).
  # If b is a non-trivial blossom,
  # dualvar[b] = z(b) where z(b) is b's variable in the dual optimization
  # problem.
  dualvar <- rep(c(maxweight, 0), each = nvertex)

  # If allowedge[k] is true, edge k has zero slack in the optimization
  # problem; if allowedge[k] is false, the edge's slack may or may not
  # be zero.
  allowedge <- rep(FALSE, nedge)

  # Queue of newly discovered S-vertices.
  queue <- c()

  # Return 2 * slack of edge k (does not work inside blossoms).
  slack <- function(idx) {
    i <- edges[2 * idx - 1]
    j <- edges[2 * idx]
    w <- weights[idx]
    return(dualvar[i] + dualvar[j] - 2 * w)
  }

  # Generate the leaf vertices of a blossom.
  blossomLeaves <- function(b) {
    rtn <- c()
    if (1 <= b && b <= nvertex) {
      rtn[length(rtn) + 1] <- b
    } else {
      for (t in blossomchilds[[b]]) {
        if (1 <= t && t <= nvertex) {
          rtn[length(rtn) + 1] <- t
        } else {
          tmp <- blossomLeaves(t)
          for (v in tmp) {
            rtn[length(rtn) + 1] <- v
          }
        }
      }
    }
    return(rtn)
  }

  # Assign label tl to the top-level blossom containing vertex w
  # and record the fact that w was reached through the edge with
  # remote endpoint p.
  assignLabel <- function(w, tl, p) {
    if (DEB) {
      DEBUG(sprintf("assignLabel(%d, %d, %d)", w, tl, p))
    }
    b <- inblossom[w]
    stopifnot(label[w] == 0 && label[b] == 0)
    label[w] <<- tl
    label[b] <<- tl
    labelend[w] <<- p
    labelend[b] <<- p
    bestedge[w] <<- -1
    bestedge[b] <<- -1
    if (tl == 1) {
      # b became an S-vertex/blossom; add it(s vertices) to the queue.
      for (v in blossomLeaves(b)) {
        queue[length(queue) + 1] <<- v
      }
      if (DEB) {
        DEBUG("PUSH ")
        DEBUG(blossomLeaves(b))
      }
    } else if (tl == 2) {
      # b became a T-vertex/blossom; assign label S to its mate.
      # (If b is a non-trivial blossom, its base is the only vertex
      # with an external mate.)
      base <- blossombase[b]
      stopifnot(mate[base] > 0)
      assignLabel(endpoint[mate[base]], 1, bitwXor(mate[base] - 1, 1) + 1)
    }
  }

  # Track back from vertices v and w to discover either a new blossom
  # or an augmenting path. Return the base vertex of the new blossom or -1.
  scanBlossom <- function(v, w) {
    if (DEB) {
      DEBUG(sprintf("scanBlossom(%d, %d)", v, w))
    }
    # Trace back from v and w, placing breadcrumbs as we go.
    mypath <<- c()
    base <- -1
    while (v != -1 || w != -1) {
      # Look for a breadcrumb in v's blossom or put a new breadcrumb.
      b <- inblossom[v]
      if (label[b] %% 8 >= 4) {
        base <- blossombase[b]
        break
      }
      stopifnot(label[b] == 1)
      mypath[length(mypath) + 1] <<- b
      label[b] <<- 5
      # Trace one step back
      stopifnot(labelend[b] == mate[blossombase[b]])
      if (labelend[b] == -1) {
        # The base of blossom b is single; stop tracing this path.
        v <- -1
      } else {
        v <- endpoint[labelend[b]]
        b <- inblossom[v]
        stopifnot(label[b] == 2)
        # b is a T-blossom; trace one more step back.
        stopifnot(labelend[b] > 0)
        v <- endpoint[labelend[b]]
      }
      # Swap v and w so that we alternate between both paths.
      if (w != -1) {
        tmp <- v
        v <- w
        w <- tmp
      }
    }
    # Remove breadcrumbs.
    for (b in mypath) {
      label[b] <<- 1
    }
    # Return base vertex, if we found one.
    return(base)
  }

  # Construct a new blossom with given base, containing edge k which
  # connects a pair of S vertices. Label the new blossom as S; set its dual
  # variable to zero; relabel its T-vertices to S and add them to the queue.
  addBlossom <- function(base, k) {
    v <- edges[2 * k - 1]
    w <- edges[2 * k]
    wt <- weights[k]
    bb <- inblossom[base]
    bv <- inblossom[v]
    bw <- inblossom[w]
    # Create blossom.
    b <- unusedblossoms[length(unusedblossoms)]
    # TODO(bowendeng): Make this more efficient
    if (length(unusedblossoms) == 1) {
      unusedblossoms <<- c()
    } else {
      unusedblossoms <<- unusedblossoms[1:(length(unusedblossoms) - 1)]
    }
    if (DEB) {
      DEBUG(sprintf("addBlossom(%d, %d) (v=%d w=%d) -> %d", base, k, v, w, b))
    }
    blossombase[b] <<- base
    blossomparent[b] <<- -1
    blossomparent[bb] <<- b
    # Make list of sub-blossoms and their interconnecting edge endpoints.
    blossomchilds[[b]] <<- list()
    mypath <- c()
    blossomendps[[b]] <<- list()
    endps <- c()
    # Trace back from v to base.
    while (bv != bb) {
      # Add bv to the new blossom.
      blossomparent[bv] <<- b
      mypath[length(mypath) + 1] <- bv
      endps[length(endps) + 1] <- labelend[bv]
      stopifnot(label[bv] == 2 ||
        (label[bv] == 1 && labelend[bv] == mate[blossombase[bv]]))
      # Trace one step back.
      stopifnot(labelend[bv] > 0)
      v <- endpoint[labelend[bv]]
      bv <- inblossom[v]
    }
    # Reverse lists, add endpoint that connects the pair of S vertices.
    mypath[length(mypath) + 1] <- bb
    mypath <- rev(mypath)
    blossomchilds[[b]] <<- mypath
    endps <- rev(endps)
    endps[length(endps) + 1] <- 2 * k - 1
    blossomendps[[b]] <- endps
    # Trace back from w to base.
    while (bw != bb) {
      # Add bw to the new blossom.
      blossomparent[bw] <<- b
      mypath[length(mypath) + 1] <- bw
      blossomchilds[[b]] <<- mypath
      endps[length(endps) + 1] <- bitwXor(labelend[bw] - 1, 1) + 1
      blossomendps[[b]] <- endps
      stopifnot(label[bw] == 2 ||
        (label[bw] == 1 && labelend[bw] == mate[blossombase[bw]]))
      # Trace one step back.
      stopifnot(labelend[bw] > 0)
      w <- endpoint[labelend[bw]]
      bw <- inblossom[w]
    }
    # Set label to S.
    stopifnot(label[bb] == 1)
    label[b] <<- 1
    labelend[b] <<- labelend[bb]
    # Set dual variable to zero.
    dualvar[b] <<- 0
    # Relabel vertices.
    for (v in blossomLeaves(b)) {
      if (label[inblossom[v]] == 2) {
        # This T-vertex now turns into an S-vertex because it becomes
        # part of an S-blossom; add it to the queue.
        queue[length(queue) + 1] <<- v
      }
      inblossom[v] <<- b
    }
    # Compute blossombestedges[[b]]
    bestedgeto <<- rep(-1, 2 * nvertex)
    for (bv in mypath) {
      if (length(blossombestedges[[bv]]) == 0) {
        # This subblossom does not have a list of least-slack edges;
        # get the information from the vertices.
        nblists <- lapply(blossomLeaves(bv), function(v) {
          sapply(neighbend[[v]], function(i) floor((i + 1) / 2))
        })
      } else {
        # Walk this subblossom's least-slack edges.
        nblists <- list()
        nblists[[1]] <- blossombestedges[[bv]]
      }
      for (idx in 1:length(nblists)) {
        nblist <- nblists[[idx]]
        if (length(nblist) == 0 || is.na(nblist)) {
          nblist <- list()
        }
        for (k in nblist) {
          i <- edges[2 * k - 1]
          j <- edges[2 * k]
          wt <- weights[k]
          if (inblossom[j] == b) {
            tmp <- i
            i <- j
            j <- tmp
          }
          bj <- inblossom[j]
          if (bj != b && label[bj] == 1 &&
            (bestedgeto[bj] == -1 ||
              slack(k) < slack(bestedgeto[bj]))) {
            bestedgeto[bj] <<- k
          }
        }
        # Forget about least-slack edges of the subblossom.
        blossombestedges[[bv]] <<- list()
        bestedge[bv] <<- -1
      }
      blossombestedges[[b]] <<- list()
      for (k in bestedgeto) {
        if (k != -1) {
          blossombestedges[[b]][length(blossombestedges[[b]]) + 1] <<- k
        }
      }
      # Select bestedge[b].
      bestedge[b] <<- -1
      for (k in blossombestedges[[b]]) {
        if (bestedge[b] == -1 || slack(k) < slack(bestedge[b])) {
          bestedge[b] <<- k
        }
      }
      if (DEB) {
        DEBUG(paste0(sprintf("blossomchilds[%d]=", b), sprintf("%s", blossomchilds[[b]]), collapse = ", "))
      }
    }
    blossomendps <<- blossomendps
  }

  # Expand the given top-level blossom.
  expandBlossom <- function(b, endstage) {
    if (DEB) {
      DEBUG(paste0(sprintf("expandBlossom(%d, %d) ", b, endstage), sprintf("%s", blossomchilds[[b]]), collapse = ", "))
    }
    # Convert sub-blossoms into top-level blossoms.
    for (s in blossomchilds[[b]]) {
      blossomparent[s] <<- -1
      if (s <= nvertex) {
        inblossom[s] <<- s
      } else if (endstage && dualvar[s] == 0) {
        # Recursively expand this sub-blossom.
        expandBlossom(s, endstage)
      } else {
        for (v in blossomLeaves(s)) {
          inblossom[v] <<- s
        }
      }
    }
    # If we expand a T-blossom during a stage, its sub-blossoms must be
    # relabeled.
    if ((!endstage) && label[b] == 2) {
      # Start at the sub-blossom through which the expanding
      # blossom obtained its label, and relabel sub-blossoms until
      # we reach the base.
      # Figure out through which sub-blossom the expanding blossom
      # obtained its label initially.
      stopifnot(labelend[b] > 0)
      entrychild <- inblossom[endpoint[bitwXor(labelend[b] - 1, 1) + 1]]
      # Decide in which direction we will go round the blossom.
      j <- which(blossomchilds[[b]] == entrychild)[1] - 1
      if (j %% 2 == 1) {
        # Start index is odd; go forward and wrap.
        j <- j - length(blossomchilds[[b]])
        jstep <- 1
        endptrick <- 0
      } else {
        # Start index is even; go backward.
        jstep <- -1
        endptrick <- 1
      }
      # Move along the blossom until we get to the base.
      p <- labelend[b]
      while (j != 0) {
        # Relabel the T-sub-blossom.
        label[endpoint[bitwXor(p - 1, 1) + 1]] <<- 0
        if (j - endptrick >= 0) {
          tmp <- blossomendps[[b]][[j + 1 - endptrick]]
        } else {
          tmp <- blossomendps[[b]][[length(blossomendps[[b]]) + 1 + j - endptrick]]
        }
        label[endpoint[bitwXor(bitwXor(tmp - 1, endptrick), 1) + 1]] <<- 0
        assignLabel(endpoint[bitwXor(p - 1, 1) + 1], 2, p)
        # Step to the next S-sub-blossom and note its forward endpoint.
        if (j - endptrick >= 0) {
          tmp <- blossomendps[[b]][[j + 1 - endptrick]]
        } else {
          tmp <- blossomendps[[b]][[length(blossomendps[[b]]) + 1 + j - endptrick]]
        }
        allowedge[floor((tmp + 1) / 2)] <<- TRUE
        j <- j + jstep
        if (j - endptrick >= 0) {
          tmp <- blossomendps[[b]][[j + 1 - endptrick]]
        } else {
          tmp <- blossomendps[[b]][[length(blossomendps[[b]]) + 1 + j - endptrick]]
        }
        p <- bitwXor(tmp - 1, endptrick) + 1
        # Step to the next T-sub-blossom.
        allowedge[floor((p + 1) / 2)] <<- TRUE
        j <- j + jstep
      }
      # Relabel the base T-sub-blossom WITHOUT stepping through to
      # its mate (so don't call assignLabel).
      if (j >= 0) {
        bv <- blossomchilds[[b]][[j + 1]]
      } else {
        bv <- blossomchilds[[b]][[length(blossomchilds[[b]]) + 1 + j]]
      }
      label[endpoint[bitwXor(p - 1, 1) + 1]] <<- 2
      label[bv] <<- 2
      labelend[endpoint[bitwXor(p - 1, 1) + 1]] <<- p
      labelend[bv] <<- p
      bestedge[bv] <<- -1
      # Continue along the blossom until we get back to entrychild.
      j <- j + jstep
      while (TRUE) {
        # Examine the vertices of the sub-blossom to see whether
        # it is reachable from a neighboring S-vertex outside the expanding blossom.
        if (j >= 0) {
          tmp <- blossomchilds[[b]][[j + 1]]
        } else {
          tmp <- blossomchilds[[b]][[length(blossomchilds[[b]]) + 1 + j]]
        }
        if (tmp == entrychild) {
          break
        }
        bv <- tmp
        if (label[bv] == 1) {
          # This sub-blossom just got label S through one of its
          # neighbors; leave it.
          j <- j + jstep
          next
        }
        for (v in blossomLeaves(bv)) {
          if (label[v] != 0) {
            break
          }
        }
        # If the sub-blossom contains a reachable vertex, assign
        # label T to the sub-blossom.
        if (label[v] != 0) {
          stopifnot(label[v] == 2)
          stopifnot(inblossom[v] == bv)
          label[v] <<- 0
          label[endpoint[mate[blossombase[bv]]]] <<- 0
          assignLabel(v, 2, labelend[v])
        }
        j <- j + jstep
      }
    }
    # Recycle the blossom number.
    label[b] <<- -1
    labelend[b] <<- -1
    blossomchilds[[b]] <<- list()
    blossomendps[[b]] <<- list()
    blossombase[b] <<- -1
    blossombestedges[[b]] <<- list()
    bestedge[b] <<- -1
    unusedblossoms[length(unusedblossoms) + 1] <<- b
  }

  # Swap matched/unmatched edges over an alternating path through blossom b
  # between vertex v and the base vertex. Keep blossom bookkeeping consistent.
  augmentBlossom <- function(b, v) {
    if (DEB) {
      DEBUG(sprintf("augmentBlossom(%d, %d)", b, v))
    }
    # Bubble up through the blossom tree from vertex v to an immediate
    # sub-blossom of b.
    t <- v
    while (blossomparent[t] != b) {
      t <- blossomparent[t]
    }
    # Recursively deal with the first sub-blossom.
    if (t > nvertex) {
      augmentBlossom(t, v)
    }
    # Decide in which direction we will go round the blossom.
    i <- j <- which(blossomchilds[[b]] == t)[1] - 1
    if (i %% 2 == 1) {
      # Start index is odd; go forward and wrap.
      j <- j - length(blossomchilds[[b]])
      jstep <- 1
      endptrick <- 0
    } else {
      # Start index is even; go backward.
      jstep <- -1
      endptrick <- 1
    }
    # Move along the blossom until we get to the base.
    while (j != 0) {
      # Step to the next sub-blossom and augment it recursively.
      j <- j + jstep
      if (j >= 0) {
        t <- blossomchilds[[b]][[j + 1]]
      } else {
        t <- blossomchilds[[b]][length(blossomchilds[[b]]) + 1 + j]
      }
      if (j - endptrick >= 0) {
        p <- bitwXor(blossomendps[[b]][[j + 1 - endptrick]] - 1, endptrick) + 1
      } else {
        p <- bitwXor(blossomendps[[b]][[length(blossomendps[[b]]) + 1 + j - endptrick]] - 1, endptrick) + 1
      }
      if (t > nvertex) {
        augmentBlossom(t, endpoint[p])
      }
      # Step to the next sub-blossom and augment it recursively.
      j <- j + jstep
      if (j >= 0) {
        t <- blossomchilds[[b]][[j + 1]]
      } else {
        t <- blossomchilds[[b]][[length(blossomchilds[[b]]) + 1 + j]]
      }
      if (t > nvertex) {
        augmentBlossom(t, endpoint[bitwXor(p - 1, 1) + 1])
      }
      # Match the edge connecting those sub-blossoms.
      mate[endpoint[p]] <<- bitwXor(p - 1, 1) + 1
      mate[endpoint[bitwXor(p - 1, 1) + 1]] <<- p
      if (DEB) {
        DEBUG(sprintf("PAIR %d %d (k=%d)", endpoint[p], endpoint[bitwXor(p - 1, 1) + 1], floor((p + 1) / 2)))
      }
    }
    # Rotate the list of sub-blossoms to put the new base at the front.
    orig_list <- blossomchilds[[b]]
    new_list <- list()
    for (idx in 1:length(orig_list)) {
      new_list[[(idx - i - 1) %% length(orig_list) + 1]] <- orig_list[[idx]]
    }
    blossomchilds[[b]] <<- new_list
    orig_list <- blossomendps[[b]]
    new_list <- list()
    for (idx in 1:length(orig_list)) {
      new_list[[(idx - i - 1) %% length(orig_list) + 1]] <- orig_list[[idx]]
    }
    blossomendps[[b]] <<- new_list
    blossombase[b] <<- blossombase[blossomchilds[[b]][[1]]]
    stopifnot(blossombase[b] == v)
  }

  # Swap matched/unmatched edges over an alternating path between two
  # single vertices. The augmenting path runs through edge k, which
  # connects a pair of S vertices.
  augmentMatching <- function(k) {
    v <- edges[2 * k - 1]
    w <- edges[2 * k]
    wt <- weights[k]
    if (DEB) {
      DEBUG(sprintf("augmentMatching(%d) (v=%d w=%d)", k, v, w))
      DEBUG(sprintf("PAIR %d %d (k=%d)", v, w, k))
    }
    for (idx in 1:2) {
      s <- c(v, w)[idx]
      p <- c(2 * k, 2 * k - 1)[idx]
      # Match vertex s to remote endpoint p. Then trace back from s
      # until we find a single vertex, swapping matched and unmatched
      # edges as we go.
      while (TRUE) {
        bs <- inblossom[s]
        stopifnot(label[bs] == 1)
        stopifnot(labelend[bs] == mate[blossombase[bs]])
        # Augment through the S-blossom from s to base.
        if (bs > nvertex) {
          augmentBlossom(bs, s)
        }
        # Update mate[s]
        mate[s] <<- p
        # Trace one step back.
        if (labelend[bs] == -1) {
          # Reached single vertex; stop.
          break
        }
        t <- endpoint[labelend[bs]]
        bt <- inblossom[t]
        stopifnot(label[bt] == 2)
        # Trace one step back.
        stopifnot(labelend[bt] > 0)
        s <- endpoint[labelend[bt]]
        j <- endpoint[bitwXor(labelend[bt] - 1, 1) + 1]
        # Augment through the T-blossom from j to base.
        stopifnot(blossombase[bt] == t)
        if (bt > nvertex) {
          augmentBlossom(bt, j)
        }
        # Update mate[j]
        mate[j] <<- labelend[bt]
        # Keep the opposite endpoint;
        # it will be assigned to mate[s] in the next step.
        p <- bitwXor(labelend[bt] - 1, 1) + 1
        if (DEB) {
          DEBUG(sprintf("PAIR %d %d (k=%d)", s, t, floor((p + 1) / 2)))
        }
      }
    }
  }

  # Verify that the optimum solution has been reached.
  verifyOptimum <- function() {
    if (maxcardinality) {
      # Vertices may have negative dual;
      # find a constant non-negative number to add to all vertex duals.
      vdualoffset <- max(0, -min(dualvar[1:nvertex]))
    } else {
      vdualoffset <- 0
    }
    # 0. all dual variables are non-negative
    stopifnot(min(dualvar[1:nvertex]) + vdualoffset >= 0)
    stopifnot(min(dualvar[(nvertex + 1):(2 * nvertex)]) >= 0)
    # 0. all edges have non-negative slack and
    # 1. all matched edges have zero slack;
    for (k in 1:nedge) {
      i <- edges[2 * k - 1]
      j <- edges[2 * k]
      wt <- weights[k]
      s <- dualvar[i] + dualvar[j] - 2 * wt
      iblossoms <- c(i)
      jblossoms <- c(j)
      while (blossomparent[iblossoms[length(iblossoms)]] != -1) {
        iblossoms[length(iblossoms) + 1] <- blossomparent[iblossoms[length(iblossoms)]]
      }
      while (blossomparent[jblossoms[length(jblossoms)]] != -1) {
        jblossoms[length(jblossoms) + 1] <- blossomparent[jblossoms[length(jblossoms)]]
      }
      iblossoms <- rev(iblossoms)
      jblossoms <- rev(jblossoms)
      for (idx in 1:length(iblossoms)) {
        bi <- iblossoms[idx]
        bj <- jblossoms[idx]
        if (bi != bj) {
          break
        }
        s <- s + 2 * dualvar[bi]
      }
      # stopifnot(s >= 0)
      if (floor((mate[i] + 1) / 2) == k || floor((mate[j] + 1) / 2) == k) {
        stopifnot(floor((mate[i] + 1) / 2) == k && floor((mate[j] + 1) / 2) == k)
        # stopifnot(s >= 0)
      }
    }
    # 2. all single vertices have zero dual value;
    for (v in 1:nvertex) {
      stopifnot(mate[v] >= 0 || dualvar[v] + vdualoffset == 0)
    }
    # 3. all blossoms with positive dual value are full.
    for (b in (nvertex + 1):(2 * nvertex)) {
      if (blossombase[b] >= 0 && dualvar[b] > 0) {
        stopifnot(length(blossomendps[[b]]) %% 2 == 1)
        for (idx in 1:length(blossomendps[[b]])) {
          if (idx %% 2 == 0) {
            p <- blossomendps[[b]][[idx]]
            stopifnot(mate[endpoint[p]] == bitwXor(p - 1, 1) + 1)
            stopifnot(mate[endpoint[bitwXor(p - 1, 1) + 1]] == p)
          }
        }
      }
    }
    # Ok.
  }

  # Check optimized delta2 against a trivial computation.
  checkDelta2 <- function() {
    for (v in 1:nvertex) {
      if (label[inblossom[v]] == 0) {
        bd <- NA
        bk <- -1
        for (p in neighbend[[v]]) {
          k <- floor((p + 1) / 2)
          w <- endpoint[p]
          if (label[inblossom[w]] == 1) {
            d <- slack(k)
            if (bk == -1 || d < bd) {
              bk <- k
              bd <- d
            }
          }
        }
        if (DEB && (bestedge[v] != -1 || bk != -1) && (bestedge[v] == -1 || bd != slack(bestedge[v]))) {
          DEBUG(paste0("v=", v, " bk=", bk, " bd=", bd, " bestedge=", bestedge[v], " slack=", slack(bestedge[v])))
        }
        stopifnot((bk == -1 && bestedge[v] == -1) || (bestedge[v] != -1 && bd == slack(bestedge[v])))
      }
    }
  }

  # Check optimized delta3 against a trivial computation.
  checkDelta3 <- function() {
    bk <- -1
    bd <- NA
    tbk <- -1
    tbd <- NA
    for (b in 1:(2 * nvertex)) {
      if (blossomparent[b] == -1 && label[b] == 1) {
        for (v in blossomLeaves(b)) {
          for (p in neighbend[[v]]) {
            k <- floor((p + 1) / 2)
            w <- endpoint[p]
            if (inblossom[w] != b && label[inblossom[w]] == 1) {
              d <- slack(k)
              if (bk == -1 || d < bd) {
                bk <- k
                bd <- d
              }
            }
          }
        }
        if (bestedge[b] != -1) {
          i <- edges[2 * bestedge[b] - 1]
          j <- edges[2 * bestedge[b]]
          wt <- weights[bestedge[b]]
          stopifnot(inblossom[i] == b || inblossom[j] == b)
          stopifnot(inblossom[i] != b || inblossom[j] != b)
          stopifnot(label[inblossom[i]] == 1 && label[inblossom[j]] == 1)
          if (tbk == -1 || slack(bestedge[b]) < tbd) {
            tbk <- bestedge[b]
            tbd <- slack(bestedge[b])
          }
        }
      }
    }
    if (DEB && ((is.na(bd) && !is.na(tbd)) ||
      (!is.na(bd) && is.na(tbd)) ||
      (!is.na(bd) && !is.na(tbd) && bd != tbd))) {
      DEBUG(sprintf("bk=%d tbk=%d bd=%s tbd=%s", bk, tbk, bd, tbd))
    }
    stopifnot((is.na(bd) && is.na(tbd)) || (!is.na(bd) && !is.na(tbd) && bd == tbd))
  }

  # Main loop: continue until no further improvement is possible.
  for (t in 1:nvertex) {
    # Each iteration of this loop is a "stage".
    # A stage finds an augmenting path and uses that to improve
    # the matching.
    if (DEB) {
      DEBUG(sprintf("STAGE %d", t))
    }

    # Remove labels from top-level blossoms/vertices.
    label <- rep(0, 2 * nvertex)

    # Forget all about least-slack edges.
    bestedge <- rep(-1, 2 * nvertex)
    for (idx in (nvertex + 1):(2 * nvertex)) {
      blossombestedges[[idx]] <- list()
    }

    # Loss of labeling means that we can not be sure that currently
    # allowable edges remain allowable throughout this stage.
    allowedge <- rep(FALSE, nedge)

    # Make queue empty.
    queue <- c()

    # Label single blossoms/vertices with S and put them in the queue.
    for (v in 1:nvertex) {
      if (mate[v] == -1 && label[inblossom[v]] == 0) {
        assignLabel(v, 1, -1)
      }
    }

    # Loop until we succeed in augmenting the matching.
    augmented <- 0
    while (TRUE) {
      # Each iteration of this loop is a "substage".
      # A substage tries to find an augmenting path;
      # if found, the path is used to improve the matching and
      # the stage ends. If there is no augmenting path, the
      # primal-dual method is used to pump some slack out of
      # the dual variables.
      if (DEB) {
        DEBUG("SUBSTAGE")
      }

      # Continue labeling until all vertices which are reachable
      # through an alternating path have got a label.
      while (length(queue) > 0 && !augmented) {
        # Take an S vertex from the queue.
        v <- queue[length(queue)]
        if (length(queue) == 1) {
          queue <- c()
        } else {
          queue <- queue[1:(length(queue) - 1)]
        }
        if (DEB) {
          DEBUG(sprintf("POP v=%d", v))
        }
        stopifnot(label[inblossom[v]] == 1)

        # Scan its neighbors:
        for (p in neighbend[[v]]) {
          k <- floor((p + 1) / 2)
          w <- endpoint[p]
          # w is a neighbor to v
          if (inblossom[v] == inblossom[w]) {
            # this edge is internal to a blossom; ignore it
            next
          }
          if (!allowedge[k]) {
            kslack <- slack(k)
            if (kslack <= 0) {
              # edge k has zero slack => it is allowable
              allowedge[k] <- TRUE
            }
          }
          if (allowedge[k]) {
            if (label[inblossom[w]] == 0) {
              # (C1) w is a free vertex;
              # label w with T and label its mate with S (R12).
              assignLabel(w, 2, bitwXor(p - 1, 1) + 1)
            } else if (label[inblossom[w]] == 1) {
              # (C2) w is an S-vertex (not in the same blossom);
              # follow back-links to discover either an
              # augmenting path or a new blossom.
              base <- scanBlossom(v, w)
              if (base >= 0) {
                # Found a new blossom; add it to the blossom
                # bookkeeping and turn it into an S-blossom.
                addBlossom(base, k)
              } else {
                # Found an augmenting path; augment the
                # matching and end this stage.
                augmentMatching(k)
                augmented <- 1
                break
              }
            } else if (label[w] == 0) {
              # w is inside a T-blossom, but w itself has not
              # yet been reached from outside the blossom;
              # mark it as reached (we need this to relabel
              # during T-blossom expansion).
              stopifnot(label[inblossom[w]] == 2)
              label[w] <- 2
              labelend[w] <- bitwXor(p - 1, 1) + 1
            }
          } else if (label[inblossom[w]] == 1) {
            # keep track of the least-slack non-allowable edge to
            # a different S-blossom.
            b <- inblossom[v]
            if (bestedge[b] == -1 || kslack < slack(bestedge[b])) {
              bestedge[b] <- k
            }
          } else if (label[w] == 0) {
            # w is a free vertex (or an unreached vertex inside
            # a T-blossom) but we can not reach it yet;
            # keep track of the least-slack edge that reaches w.
            if (bestedge[w] == -1 || kslack < slack(bestedge[w])) {
              bestedge[w] <- k
            }
          }
        }
      }
      if (augmented) {
        break
      }
      # There is no augmenting path under these constraints;
      # compute delta and reduce slack in the optimization problem.
      # (Note that our vertex dual variables, edge slacks and delta's
      # are pre-multiplied by two.)
      deltatype <- -1
      delta <- deltaedge <- deltablossom <- c()

      # Verify data structures for delta2/delta3 computation.
      if (CHECK_DELTA) {
        checkDelta2()
        checkDelta3()
      }

      # Compute delta1: the minimum value of any vertex dual.
      if (!maxcardinality) {
        deltatype <- 1
        delta <- min(dualvar[1:nvertex])
      }

      # Compute delta2: the minimum slack on any edge between
      # an S-vertex and a free vertex.
      for (v in 1:nvertex) {
        if (label[inblossom[v]] == 0 && bestedge[v] != -1) {
          d <- slack(bestedge[v])
          if (deltatype == -1 || d < delta) {
            delta <- d
            deltatype <- 2
            deltaedge <- bestedge[v]
          }
        }
      }

      # Compute delta3: half the minimum slack on any edge between
      # a pair of S-blossoms.
      for (b in 1:(2 * nvertex)) {
        if (blossomparent[b] == -1 && label[b] == 1 && bestedge[b] != -1) {
          kslack <- slack(bestedge[b])
          if (floor(kslack) == kslack) {
            stopifnot(kslack %% 2 == 0)
            d <- floor(kslack / 2)
          } else {
            d <- kslack / 2
          }
          if (deltatype == -1 || d < delta) {
            delta <- d
            deltatype <- 3
            deltaedge <- bestedge[b]
          }
        }
      }

      # Compute delta4: minimum z variable of any T-blossom.
      for (b in (nvertex + 1):(2 * nvertex)) {
        if (blossombase[b] >= 0 && blossomparent[b] == -1 &&
          label[b] == 2 &&
          (deltatype == -1 || dualvar[b] < delta)) {
          delta <- dualvar[b]
          deltatype <- 4
          deltablossom <- b
        }
      }

      if (deltatype == -1) {
        # No further improvement possible; max-cardinality optimum
        # reached. Do a final delta update to make the optimum
        # verifyable.
        stopifnot(maxcardinality)
        deltatype <- 1
        delta <- max(0, min(dualvar[1:nvertex]))
      }

      # Update dual variables according to delta.
      for (v in 1:nvertex) {
        if (label[inblossom[v]] == 1) {
          # S-vertex: 2*u = 2*u - 2*delta
          dualvar[v] <- dualvar[v] - delta
        } else if (label[inblossom[v]] == 2) {
          # T-vertex: 2*u = 2*u + 2*delta
          dualvar[v] <- dualvar[v] + delta
        }
      }
      for (b in (nvertex + 1):(2 * nvertex)) {
        if (blossombase[b] >= 0 && blossomparent[b] == -1) {
          if (label[b] == 1) {
            # top-level S-blossom: z = z + 2*delta
            dualvar[b] <- dualvar[b] + delta
          } else if (label[b] == 2) {
            # top-level T-blossom: z = z - 2*delta
            dualvar[b] <- dualvar[b] - delta
          }
        }
      }

      # Take action at the point where minimum delta occurred.
      if (DEB) {
        DEBUG(sprintf("delta%d=%f", deltatype, delta))
      }
      if (deltatype == 1) {
        # No further improvement possible; optimum reached.
        break
      } else if (deltatype == 2) {
        # Use the least-slack edge to continue the search.
        allowedge[deltaedge] <- TRUE
        i <- edges[2 * deltaedge - 1]
        j <- edges[2 * deltaedge]
        wt <- weights[deltaedge]
        if (label[inblossom[i]] == 0) {
          tmp <- i
          i <- j
          j <- tmp
        }
        stopifnot(label[inblossom[i]] == 1)
        queue[length(queue) + 1] <- i
      } else if (deltatype == 3) {
        # Use the least-slack edge to continue the search.
        allowedge[deltaedge] <- TRUE
        i <- edges[2 * deltaedge - 1]
        j <- edges[2 * deltaedge]
        wt <- weights[deltaedge]
        stopifnot(label[inblossom[i]] == 1)
        queue[length(queue) + 1] <- i
      } else if (deltatype == 4) {
        # Expand the least-z blossom.
        expandBlossom(deltablossom, FALSE)
      }

      # End of a this substage.
    }

    # Stop when no more augmenting path can be found.
    if (!augmented) {
      break
    }

    # End of a stage; expand all S-blossoms which have dualvar = 0.
    for (b in (nvertex + 1):(2 * nvertex)) {
      if (blossomparent[b] == -1 && blossombase[b] >= 0 &&
        label[b] == 1 && dualvar[b] == 0) {
        expandBlossom(b, TRUE)
      }
    }
  }

  # Verify that we reached the optimum solution.
  if (CHECK_OPTIMUM) {
    verifyOptimum()
  }

  # Transform mate[] such that mate[v] is the vertex to which v is paired.
  for (v in 1:nvertex) {
    if (mate[v] >= 0) {
      mate[v] <- endpoint[mate[v]]
    }
  }
  for (v in 1:nvertex) {
    stopifnot(mate[v] == -1 || mate[mate[v]] == v)
  }

  matching_size <- floor(sum(mate != -1) / 2)
  matching <- rep(NA, matching_size * 2)
  matching_weight <- rep(NA, matching_size)
  # matched[i] is TRUE when is is matched
  matched <- rep(FALSE, nvertex)
  matching_num <- 1
  for (v in 1:nvertex) {
    if (!matched[v] && mate[v] != -1) {
      w <- mate[v]
      # v must be smaller than w
      stopifnot(v < w)
      wt <- edge_weight[v, w]
      matching[matching_num * 2 - 1] <- v
      matching[matching_num * 2] <- w
      matching_weight[matching_num] <- wt
      matched[v] <- matched[w] <- TRUE
      matching_num <- matching_num + 1
    }
  }
  return(list(matching_size = matching_size, matching_weight = matching_weight, matching = matching, mate = mate))
}
