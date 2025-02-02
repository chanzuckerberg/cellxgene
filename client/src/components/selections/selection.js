import React from "react";
import cls from "./selection.css";


export class Selection extends React.PureComponent {
  render() {
    const { name, indexes, onMouseOver, onMouseLeave } = this.props;

    return (
      <div
        className={cls.root}
        onFocus={onMouseOver}
        onMouseOver={onMouseOver}
        onMouseLeave={onMouseLeave}
      >
        <span style={{ fontWeight: "bold" }}>{name}</span>
        <span>{indexes.length} cells</span>
      </div>
    );
  }
}
