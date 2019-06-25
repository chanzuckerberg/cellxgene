// jshint esversion: 6
import React from "react";

import styles from "./container.css";

const Container = props => {
  const { children } = props;
  return <div className={styles.container}>{children}</div>;
};

export default Container;
