import React, { PureComponent } from "react";
import { Drawer } from "@blueprintjs/core";

export default class InfoDrawer extends PureComponent {
  constructor(props) {
    super(props);

    this.state = {
      isOpen: false,
    };
  }

  handleClick = () => {
    this.setState((state) => {
      return { isOpen: !state.isOpen };
    });
  };

  handleKeyPress = (e) => {
    if (e.key === "Enter") {
      this.handleCategoryClick();
    }
  };

  render() {
    const { isOpen } = this.state;
    const {
      position,
      aboutURL,
      datasetTitle,
      title,
      children,
      metadataField,
      theOneValue,
    } = this.props;

    return (
      <div
        role="menuitem"
        tabIndex="0"
        onKeyPress={this.handleKeyPress}
        onClick={this.handleClick}
      >
        {children}
        <Drawer {...{ isOpen, position, title }}>
          Here is a test
          <a
            style={{ width: 185 }}
            href={aboutURL}
            data-testid="header"
            target="_blank"
            rel="noopener noreferrer"
          >
            {datasetTitle}
          </a>
          (
          <div style={{ marginBottom: 10, marginTop: 4 }}>
            <span style={{ maxWidth: 150, fontWeight: 700 }}>
              {metadataField}
            </span>
            <span style={{ maxWidth: 150 }}>{`: ${theOneValue}`}</span>
          </div>
          );
        </Drawer>
      </div>
    );
  }
}
