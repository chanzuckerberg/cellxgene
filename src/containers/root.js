import React from 'react';

import styles from './root.css';

const Root = props => (
  <div className={styles.root}>
    {props.children}
  </div>
);

export default Root;
