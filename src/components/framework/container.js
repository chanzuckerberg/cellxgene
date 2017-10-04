import React from 'react';

import styles from './container.css';

const Container = props => (
  <div className={styles.container}>
    {props.children}
  </div>
);

export default Container;
