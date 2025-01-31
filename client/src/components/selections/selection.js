import React from "react";

export class Selection extends React.PureComponent {
  render() {
    const { name, indexes, onMouseOver, onMouseLeave } = this.props;

    return (
      <div
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          gap: "2ch",
          cursor: "pointer",
        }}
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
